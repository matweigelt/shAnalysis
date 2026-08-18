"""Audit-driven help completion. Reads help_audit findings, patches files."""
import re, sys, subprocess, os
sys.path.insert(0, os.path.join('/tmp/shx_git', 'tools'))

OPTDESC = {
 'kn': 'load Love numbers, degrees 0..nmax (user-supplied; e.g. shLowLevel.fetchLoveNumbers)',
 'hn': 'vertical-deformation Love numbers, degrees 0..nmax (user-supplied)',
 'ln': 'horizontal-deformation Love numbers, degrees 0..nmax (user-supplied)',
 'Timeout': 'network timeout [s]',
 'Proxy': 'per-call proxy URL, e.g. "http://proxy:8080" (empty: MATLAB Web Preferences)',
 'Update': 'refresh existing files (safe swap: verified before replacing)',
 'Quiet': 'suppress progress output',
 'Dest': 'target folder (empty: subfolder of shLowLevel.dataFolder)',
 'BaseURL': 'server base URL, or a local mirror folder for offline use',
 'Names': 'display labels, one per solution/series',
 'Plot': 'render the diagnostic figure(s)',
 'Mask': 'nlat x nlon logical region restriction',
 'quantity': 'physical quantity of the field values',
 'Quantity': 'physical quantity of the field values',
 'MinDegree': 'first spherical harmonic degree included',
 'Sidecar': 'write <file>.provenance.json alongside the output',
 'T0': 'reference epoch [decimal years] (NaN: mean of the epochs)',
 'idx': 'index struct from shLowLevel.shIndex matching the coefficient vectors',
 'tYears': 'epochs [decimal years] ([]: taken from the series)',
 'LatDeg': 'latitude grid [deg, geocentric]',
 'LonDeg': 'longitude grid [deg]',
 'Catalog': 'row numbers of the shLowLevel.listITSG catalogue to fetch completely',
 'Release': 'explicit ITSG release (empty: 2018 before 2017.5, operational after)',
 'List': 'pre-fetched catalogue table (avoids repeated server queries)',
 'DataFolder': 'persistent data folder, applied BEFORE any fetcher runs',
}
INDESC = {}   # filled per iteration for flagged inputs

def load():
    import importlib, help_audit
    importlib.reload(help_audit)
    return help_audit

def findings():
    r = subprocess.run(['python3', 'tools/help_audit.py'],
                       capture_output=True, text=True, cwd='/tmp/shx_git')
    return [l.strip()[2:] for l in r.stdout.splitlines() if l.startswith('  - ')]

def locate(label):
    """label -> (path, entry) using the audit parsers"""
    ha = load()
    if label.startswith('shLowLevel.'):
        p = os.path.join('/tmp/shx_git', '+shLowLevel', label.split('.')[1] + '.m')
        return p, ha.parse_file(p)
    if '.' in label and label.split('.')[0] in ('shCoefficients', 'shSeries', 'shClimatology'):
        cls, meth = label.split('.')
        p = os.path.join('/tmp/shx_git', cls + '.m')
        lines = open(p, encoding='utf-8').readlines()
        for i, L in enumerate(lines):
            if re.match(r'^\s*function\s+.*\b' + meth + r'\s*\(', L):
                e = ha.parse_from(lines, i, '')
                if e and e['name'] == meth:
                    return p, e
        raise KeyError(label)
    p = os.path.join('/tmp/shx_git', label + '.m')
    return p, ha.parse_file(p)

def helpspan(path, entry):
    lines = open(path, encoding='utf-8').readlines()
    i = entry['line'] + 1
    j = i
    while j < len(lines) and lines[j].lstrip().startswith('%'):
        j += 1
    prefix = re.match(r'^(\s*)%', lines[i]).group(1) if j > i else \
             re.match(r'^(\s*)', lines[entry['line']]).group(1)
    return lines, i, j, prefix

def insert_before_anchor(lines, i, j, prefix, block):
    """insert help lines before Example/Developed/end-of-help"""
    at = j
    for k in range(i, j):
        s = lines[k].lstrip().lstrip('%').strip()
        if re.match(r'^(Examples?\b|Developed by)', s):
            at = k
            break
    return lines[:at] + block + lines[at:], at

def fix_option_default(label, opt, default):
    path, e = locate(label)
    lines, i, j, pre = helpspan(path, e)
    tgt = None
    for k in range(i, j):
        if re.search(r'\b' + re.escape(opt) + r'\b', lines[k]):
            tgt = k; break
    dtxt = default if default else '""'
    if tgt is None:
        return fix_option_missing(label, opt, default)
    L = lines[tgt]
    L2 = re.sub(r'\b' + re.escape(opt) + r'\b', opt + ' (' + dtxt + ')', L, count=1)
    lines[tgt] = L2
    open(path, 'w', encoding='utf-8').write(''.join(lines))
    return True

def fix_option_missing(label, opt, default):
    path, e = locate(label)
    lines, i, j, pre = helpspan(path, e)
    desc = OPTDESC.get(opt, 'see arguments block')
    dtxt = default if default else '""'
    # find existing Options section to append into
    at = None
    for k in range(i, j):
        if re.match(r'^\s*%\s*Options\b', lines[k]):
            at = k + 1
            while at < j and lines[at].lstrip().startswith('%') and \
                  re.match(r'^\s*%\s{4,}', lines[at]):
                at += 1
            break
    row = pre + '%     ' + opt + ' (' + dtxt + ')  ' + desc + '\n'
    if at is not None:
        lines = lines[:at] + [row] + lines[at:]
    else:
        block = [pre + '%\n', pre + '%   Options\n', row]
        lines, _ = insert_before_anchor(lines, i, j, pre, block)
    open(path, 'w', encoding='utf-8').write(''.join(lines))
    return True

def fix_input(label, name):
    key = (label, name)
    desc = INDESC.get(key) or INDESC.get(name)
    if not desc: return False
    path, e = locate(label)
    lines, i, j, pre = helpspan(path, e)
    block = [pre + '%     ' + name + '  ' + desc + '\n']
    # append into Inputs section if any, else create before Options/Outputs
    at = None
    for k in range(i, j):
        if re.match(r'^\s*%\s*Inputs\b', lines[k]):
            at = k + 1; break
    if at is not None:
        lines = lines[:at] + block + lines[at:]
    else:
        blk = [pre + '%\n', pre + '%   Inputs\n'] + block
        for k in range(i, j):
            s = lines[k].lstrip().lstrip('%').strip()
            if re.match(r'^(Options|Outputs?)\b', s):
                lines = lines[:k] + blk + lines[k:]
                break
        else:
            lines, _ = insert_before_anchor(lines, i, j, pre, blk)
    open(path, 'w', encoding='utf-8').write(''.join(lines))
    return True

OUTDESC = {}  # (label) -> list of '%   ...' body lines (without prefix)

OUTDESC.update({
 'shLowLevel.csFromVec': ['%     cnm  (nmax+1 x nmax+1) double  cosine coefficients, C(n+1, m+1)',
                          '%     snm  (nmax+1 x nmax+1) double  sine coefficients, S(:, 1) = 0'],
 'shLowLevel.errorMap': ['%     sig     (nlat x nlon) double  1-sigma of the synthesized quantity',
                         '%     latDeg  (1, nlat) double  latitude grid used [deg]',
                         '%     lonDeg  (1, nlon) double  longitude grid used [deg]'],
 'shLowLevel.mcPropagate': ['%     out  (1,1) struct  fields: samples (K x nSamples double), mean,',
                            '%          sigma (K x 1 double)  Monte-Carlo moments of the functional'],
 'shLowLevel.opApply': ['%     Yout  (P x T) double  operator applied to every column of Y'],
 'shLowLevel.parseGraceFilename': ['%     meta  (1,1) struct  fields: center/product/release (string),',
                                   '%           epoch/year/month (double), nmax (double, NaN if absent)'],
 'shLowLevel.readSINEX': ['%     snx  (1,1) struct  fields: x (P x 1 double), idx (struct),',
                          '%          Cxx (P x P double, [] unless Only="full"), epoch (1,1 double),',
                          '%          N/b (normal equations when present)'],
 'shLowLevel.readTN14': ['%     tn  (1,1) struct  fields: t (T x 1 double, decimal years),',
                         '%         C20/C30 (T x 1 double), sigmaC20/sigmaC30 (T x 1 double)'],
 'shLowLevel.resolutionMap': ['%     res  (nlat x nlon) double  local resolution [km] of the filter'],
 'shLowLevel.shAnalysisGrid': ['%     info  (1,1) struct  fields: condest (1,1 double), nObs (1,1 double),',
                               '%           kaula (1,1 double, NaN when unregularized)'],
 'shLowLevel.shEvalGFCT': ['%     Ct  (nmax+1 x nmax+1) double  cosine coefficients at the epoch',
                           '%     St  (nmax+1 x nmax+1) double  sine coefficients at the epoch'],
 'shLowLevel.shReadGFC': ['%     model  (1,1) struct  fields: C/S/sigmaC/sigmaS (nmax+1 x nmax+1',
                          '%            double), GM/R (1,1 double), nmax (1,1 double), tide (string),',
                          '%            variable terms (trend/annual/semiannual) when present'],
 'shLowLevel.shSpectralCrossover': ['%     nCrossover    (1,1) double  first degree with err >= signal (NaN if none)',
                                    '%     degreeInterp  (1,1) double  linearly interpolated crossing degree'],
 'shLowLevel.shSynthesis': ['%     grid  (nlat x nlon) double  synthesized quantity',
                            '%     lat   (1, nlat) double  latitude grid [deg]',
                            '%     lon   (1, nlon) double  longitude grid [deg]',
                            '%     P     (1,1) struct  cached Legendre functions (pass back in for reuse)'],
 'shLowLevel.plotBasinSeries': ['%     h  (1,1) graphics handle  axes of the basin series plot'],
 'shLowLevel.plotSHMap': ['%     h  (1,1) graphics handle  axes of the map plot'],
 'runAllTests': ['%     results  (1,N) matlab.unittest.TestResult  full result array',
                 '%              (pass/fail/duration per test; log written to LogFile)'],
 'shCoefficients.crossover': ['%     n0       (1,1) double  first degree with err >= signal (NaN if none)',
                              '%     nInterp  (1,1) double  interpolated crossing degree'],
 'shClimatology.removeGIA': ['%     out  (1,1) shClimatology  copy with the GIA trend removed from the',
                             '%          trend component; history appended'],
})


def fix_outputs(label):
    path, e = locate(path_label := label) if False else locate(label)
    lines, i, j, pre = helpspan(path, e)
    cls = label.split('.')[0]
    body = OUTDESC.get(label)
    if body is None:
        outs = e['outs']
        if cls in ('shCoefficients', 'shSeries', 'shClimatology') and \
           outs in (['out'], ['obj']):
            body = ['%     ' + outs[0] + '  (1,1) ' + cls +
                    '  modified copy; the operation is appended to the']
            body += ['%              history (immutable value-class pattern)']
        else:
            return False
    block = [pre + '%\n', pre + '%   Outputs\n'] + \
            [pre + b + '\n' for b in body]
    lines, _ = insert_before_anchor(lines, i, j, pre, block)
    open(path, 'w', encoding='utf-8').write(''.join(lines))
    return True

def run():
    os.chdir('/tmp/shx_git')
    F = findings()
    done = skip = 0
    for f in F:
        m = re.match(r"(\S+): option '(\w+)' documented without its default \((.*)\)$", f)
        if m:
            done += fix_option_default(m.group(1), m.group(2),
                                       m.group(3) if m.group(3) != '?' else '')
            continue
        m = re.match(r"(\S+): option '(\w+)' undocumented$", f)
        if m:
            done += fix_option_missing(m.group(1), m.group(2), default_of(m.group(1), m.group(2)))
            continue
        m = re.match(r"(\S+): input '(\w+)' undocumented$", f)
        if m:
            ok = fix_input(m.group(1), m.group(2))
            done += ok; skip += not ok
            continue
        m = re.match(r"(\S+): no Outputs section$", f)
        if m:
            ok = fix_outputs(m.group(1))
            done += ok; skip += not ok
            continue
        skip += 1
    print('fixed:', done, 'needs curation:', skip)

def default_of(label, opt):
    _, e = locate(label)
    for a in e['args']:
        if a['name'] == opt and a['nv']:
            return a['default']
    return ''

if __name__ == '__main__':
    run()
