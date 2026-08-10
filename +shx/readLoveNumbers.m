function [LN, info] = readLoveNumbers(filename, opts)
%READLOVENUMBERS Read load Love numbers from a user-supplied table file.
%
%   LN = shx.readLoveNumbers(FILE) parses a plain-text load-Love-number
%   table (comment lines starting with %, #, ! or // are skipped; an
%   optional header line naming the columns is recognized). The toolbox
%   never hardcodes Love numbers - this reader turns YOUR loading-model
%   table (PREM, Gutenberg-Bullen, ak135, ...) into the kn/hn/ln vectors
%   the synthesis and deformation routines expect.
%
%   Recognized layouts
%     2 columns:            n  kn            (assumed; documented)
%     3+ columns w/ header: mapped by the header tokens n, k(n)('), 
%                           h(n)('), l(n)(')
%     3+ columns w/o header: AMBIGUOUS (classic Farrell tables order
%                           h l k, others k h l) - pass Columns=
%                           explicitly, e.g. Columns="n h l k".
%
%   Options
%     Columns ("auto")   explicit column meaning, tokens from {n,k,h,l},
%                        e.g. "n k h l" (overrides header detection)
%     MaxDegree ([])     verify coverage of degrees 0..MaxDegree and
%                        truncate; with sparse tables combine with Interp
%     Interp ("none")    "none" | "pchip": fill degree gaps by shape-
%                        preserving pchip in log(1+n) - Farrell-style
%                        sparse tables (n = ...,32,56,100,...) become
%                        dense. Python-validated: max abs error 1.3e-4
%                        (0.05% of curve range) reconstructing a smooth
%                        kn' curve from 15 Farrell-spaced samples.
%                        Interpolation stays inside the table span;
%                        extrapolation errors out.
%     InFrame ("")       degree-1 reference frame of the FILE, and
%     OutFrame ("")      requested output frame: "CE" | "CF" | "CM".
%                        Both or neither. Only degree 1 changes (Blewitt
%                        2003 isomorphic frames): all of h1,l1,k1 shift
%                        by the same constant d:
%                          CE->CF  d = (h1+2*l1)/3
%                          CE->CM  d = 1 + k1
%                          CF->CM  d = 1 + k1          (k1 in CF values)
%                          CM->CF  d = (h1+2*l1)/3     (in CM values)
%                        Conversions TO "CE" from CF/CM are not
%                        recoverable from the table alone and error.
%                        Validated against the published PREM values:
%                        CE (-0.290, 0.113, 0.021) -> CF (-0.269, 0.134)
%                        and the exact identity k1(CM) = -1; the
%                        CF->CM->CF roundtrip is exact.
%
%   Outputs
%     LN    struct: n     (D x 1) degrees, ascending
%                   kn    (D x 1) potential load Love numbers k'_n
%                   hn    (D x 1) vertical    (NaN if absent from file)
%                   ln    (D x 1) horizontal  (NaN if absent from file)
%     info  struct: source, columns, frame, interpolated (logical mask)
%
%   Example (G7 GNSS chain)
%     LN = shx.readLoveNumbers("prem_load.txt", Columns="n h l k", ...
%              MaxDegree=96, Interp="pchip", InFrame="CE", OutFrame="CF");
%     [up, no, ea] = g.deformation(latPts, lonPts, kn=LN.kn, hn=LN.hn, ...
%              ln=LN.ln, Mode="points");
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-07 (v2.5).
arguments
    filename {mustBeTextScalar}
    opts.Columns (1,1) string = "auto"
    opts.MaxDegree double {mustBeScalarOrEmpty} = []
    opts.Interp (1,1) string {mustBeMember(opts.Interp, ...
        ["none", "pchip"])} = "none"
    opts.InFrame (1,1) string {mustBeMember(opts.InFrame, ...
        ["", "CE", "CF", "CM"])} = ""
    opts.OutFrame (1,1) string {mustBeMember(opts.OutFrame, ...
        ["", "CE", "CF", "CM"])} = ""
end
filename = char(filename);
if ~isfile(filename)
    error('shx:readLoveNumbers:fileNotFound', 'File not found: %s', ...
        filename);
end
if xor(strlength(opts.InFrame) > 0, strlength(opts.OutFrame) > 0)
    error('shx:readLoveNumbers:badFrame', ...
        'InFrame and OutFrame must be given together (or neither).');
end

lines = readlines(filename);
vals = {}; header = "";
for k = 1:numel(lines)
    L = strtrim(lines(k));
    if L == "", continue; end
    isComment = startsWith(L, ["%", "#", "!", "//"]);
    if isComment
        % commented lines are never data, but the LAST one before the
        % first data row may carry the column header ("# n h' l' k'")
        if isempty(vals)
            header = strtrim(erase(L, ["%", "#", "!", "//"]));
        end
        continue;
    end
    toksL = split(L, whitespacePattern | "," | ";");
    toksL = toksL(strlength(toksL) > 0);
    q = str2double(toksL);
    if any(isnan(q))
        if isempty(vals), header = L; end           % bare header line
        continue;
    end
    vals{end+1} = q(:)'; %#ok<AGROW>
end
if isempty(vals)
    error('shx:readLoveNumbers:noData', ...
        'No numeric rows found in %s.', filename);
end
ncol = numel(vals{1});
if any(cellfun(@numel, vals) ~= ncol)
    error('shx:readLoveNumbers:noData', ...
        'Inconsistent column count in %s.', filename);
end
M = vertcat(vals{:});

% ------------------------------------------------------- column mapping
map = resolveColumns(opts.Columns, header, ncol, filename);
n = M(:, map.n);
if any(n < 0 | n ~= round(n))
    error('shx:readLoveNumbers:badColumns', ...
        ['Degree column contains non-integers - check the Columns= ' ...
         'mapping (got "%s").'], strjoin(map.names, ' '));
end
[n, iu] = unique(n, 'sorted');
kn = M(iu, map.k);
hn = nan(size(n)); ln = nan(size(n));
if map.h > 0, hn = M(iu, map.h); end
if map.l > 0, ln = M(iu, map.l); end

% ----------------------------------------------- coverage/interpolation
interpolated = false(size(n));
if ~isempty(opts.MaxDegree)
    N = opts.MaxDegree;
    if n(1) > 0 || n(end) < N
        error('shx:readLoveNumbers:degreeGap', ...
            ['Table spans degrees %d..%d but 0..%d requested - ' ...
             'interpolation cannot extrapolate.'], n(1), n(end), N);
    end
    want = (0:N)';
    miss = ~ismember(want, n);
    if any(miss) && opts.Interp == "none"
        error('shx:readLoveNumbers:degreeGap', ...
            ['Degrees missing from the table (first gap at n=%d). ' ...
             'Sparse Farrell-style tables need Interp="pchip".'], ...
            want(find(miss, 1)));
    end
    if any(miss)
        x = log1p(n);  xi = log1p(want);
        knI = interp1(x, kn, xi, 'pchip');
        hnI = nan(size(want)); lnI = hnI;
        if map.h > 0, hnI = interp1(x, hn, xi, 'pchip'); end
        if map.l > 0, lnI = interp1(x, ln, xi, 'pchip'); end
        interpolated = miss;
        n = want; kn = knI; hn = hnI; ln = lnI;
    else
        keep = ismember(n, want);
        n = n(keep); kn = kn(keep); hn = hn(keep); ln = ln(keep);
        interpolated = false(size(n));
    end
end

% -------------------------------------------------- degree-1 frame shift
frameNote = "as-in-file";
if strlength(opts.OutFrame) > 0 && opts.InFrame ~= opts.OutFrame
    i1 = find(n == 1, 1);
    if isempty(i1)
        error('shx:readLoveNumbers:needsDegree1', ...
            'Frame conversion needs degree 1 in the table.');
    end
    pairKey = opts.InFrame + ">" + opts.OutFrame;
    switch pairKey
        case "CE>CF"
            needHL(map, pairKey);
            d = (hn(i1) + 2 * ln(i1)) / 3;
        case "CE>CM"
            d = 1 + kn(i1);
        case "CF>CM"
            d = 1 + kn(i1);
        case "CM>CF"
            needHL(map, pairKey);
            d = (hn(i1) + 2 * ln(i1)) / 3;
        otherwise   % "CF>CE", "CM>CE"
            error('shx:readLoveNumbers:frameNotRecoverable', ...
                ['%s: the CE-frame values are not recoverable from ' ...
                 '%s-frame numbers alone (the defining shift is ' ...
                 'expressed in CE values; Blewitt 2003).'], ...
                pairKey, opts.InFrame);
    end
    kn(i1) = kn(i1) - d;
    if map.h > 0, hn(i1) = hn(i1) - d; end
    if map.l > 0, ln(i1) = ln(i1) - d; end
    frameNote = pairKey + sprintf(" (degree-1 shift %+.6g)", d);
end

LN = struct('n', n, 'kn', kn, 'hn', hn, 'ln', ln);
info = struct('source', string(filename), ...
    'columns', strjoin(map.names, ' '), ...
    'frame', frameNote, 'interpolated', interpolated);
end

function needHL(map, key)
if map.h <= 0 || map.l <= 0
    error('shx:readLoveNumbers:needsDegree1', ...
        '%s needs h and l columns in the table.', key);
end
end

function map = resolveColumns(colsOpt, header, ncol, filename)
%RESOLVECOLUMNS Map column indices for n/k/h/l from option or header.
map = struct('n', 0, 'k', 0, 'h', 0, 'l', 0, 'names', {{}});
toks = string.empty; src = "";
if colsOpt ~= "auto"
    toks = lower(strtrim(split(colsOpt, whitespacePattern | ",")));
    toks = toks(strlength(toks) > 0);
    src = "Columns option";
elseif strlength(header) > 0
    % a header is only trusted if every token cleans to n/k/h/l and the
    % count matches the data columns - prose comments are ignored
    ht = lower(regexprep(split(strtrim(header), whitespacePattern), ...
        "[_'()\d]", ""));
    ht = ht(strlength(ht) > 0);
    ok = numel(ht) == ncol && all(ismember(extractBefore(ht + "x", 2), ...
        ["n", "d", "k", "h", "l"]));      % d = "degree"
    if ok
        toks = ht; src = "header line";
    end
end
if isempty(toks)
    if ncol == 2
        toks = ["n"; "k"];
        src = "2-column default [n kn]";
    else
        error('shx:readLoveNumbers:ambiguousColumns', ...
            ['%d columns without a recognizable header: the ordering ' ...
             'is ambiguous (Farrell tables use n h l k, others ' ...
             'n k h l). Pass Columns= explicitly, e.g. ' ...
             'Columns="n h l k".'], ncol);
    end
end
if numel(toks) ~= ncol
    error('shx:readLoveNumbers:badColumns', ...
        ['Column mapping (%s: %s) does not match the %d data ' ...
         'columns of %s.'], src, strjoin(toks, ' '), ncol, filename);
end
names = cell(1, ncol);
for c = 1:ncol
    t = char(toks(c));
    switch t(1)
        case {'n', 'd'}, map.n = c; names{c} = 'n';   % n or "degree"
        case 'k', map.k = c; names{c} = 'k';
        case 'h', map.h = c; names{c} = 'h';
        case 'l', map.l = c; names{c} = 'l';
        otherwise
            names{c} = '?';                    % ignored column
    end
end
map.names = names;
if map.n == 0 || map.k == 0
    error('shx:readLoveNumbers:badColumns', ...
        'Column mapping must include n and k (got: %s).', ...
        strjoin(names, ' '));
end
end
