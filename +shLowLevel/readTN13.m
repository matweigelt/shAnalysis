function tn = readTN13(filename)
%READTN13 Read a GRACE/GRACE-FO TN-13 geocenter (degree-1) file.
%
%   TN = shLowLevel.readTN13(FILENAME) parses the official TN-13 product
%   (JPL/CSR/GFZ, e.g. TN-13_GEOC_CSR_RL06.txt): data lines
%       GRCOF2  n  m  Clm  Slm  Clm_sig  Slm_sig  begin  end
%   with n=1, m=0/1 and begin/end epochs yyyymmdd.hhmm. Line pairs
%   (m=0, m=1) sharing the same interval are combined into one monthly
%   record with the interval mid-point as epoch.
%
%   Inputs
%     filename  char/string  path (.txt, or .gz which is gunzipped)
%   Outputs
%     tn  struct:
%       epoch  (K,1) double  interval mid-points [decimal years]
%       C10, C11, S11        (K,1) double  degree-1 Stokes coefficients
%       sigC10, sigC11, sigS11 (K,1) double 1-sigma
%       t0, t1 (K,1) double  interval bounds [decimal years]
%
%   Format assumption (documented): columns are whitespace-separated in
%   the order above; header lines (anything not starting with GRCOF2)
%   are skipped. Verify against your file's yaml header on first use.
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    filename {mustBeTextScalar}
end
filename = char(filename);
if ~isfile(filename)
    error('shLowLevel:readTN13:fileNotFound', 'File not found: %s', filename);
end
if endsWith(filename, '.gz')
    tmp = tempname; mkdir(tmp);
    files = gunzip(filename, tmp);
    filename = files{1};
end
lines = readlines(filename);
rec = struct('t0', {}, 't1', {}, 'n', {}, 'm', {}, 'C', {}, 'S', {}, ...
    'sC', {}, 'sS', {});
for k = 1:numel(lines)
    L = strtrim(lines(k));
    if ~startsWith(L, "GRCOF2"), continue; end
    p = str2double(split(L));
    if numel(p) < 9 || any(isnan(p(2:9)))
        error('shLowLevel:readTN13:badLine', 'Unparseable GRCOF2 line: %s', L);
    end
    rec(end+1) = struct('t0', p(8), 't1', p(9), 'n', p(2), 'm', p(3), ...
        'C', p(4), 'S', p(5), 'sC', p(6), 'sS', p(7)); %#ok<AGROW>
end
if isempty(rec)
    error('shLowLevel:readTN13:noData', 'No GRCOF2 records found in %s.', filename);
end
key = [[rec.t0]', [rec.t1]'];
[uk, ~, iu] = unique(key, 'rows', 'stable');
K = size(uk, 1);
tn.t0 = shLowLevel.icgemDate2Year(uk(:,1));
tn.t1 = shLowLevel.icgemDate2Year(uk(:,2));
tn.epoch = (tn.t0 + tn.t1) / 2;
z = nan(K, 1);
tn.C10 = z; tn.C11 = z; tn.S11 = z;
tn.sigC10 = z; tn.sigC11 = z; tn.sigS11 = z;
for j = 1:numel(rec)
    i = iu(j);
    if rec(j).n ~= 1
        error('shLowLevel:readTN13:badDegree', ...
            'TN-13 must contain only degree-1 lines (got n=%d).', rec(j).n);
    end
    if rec(j).m == 0
        tn.C10(i) = rec(j).C;    tn.sigC10(i) = rec(j).sC;
    else
        tn.C11(i) = rec(j).C;    tn.sigC11(i) = rec(j).sC;
        tn.S11(i) = rec(j).S;    tn.sigS11(i) = rec(j).sS;
    end
end
end
