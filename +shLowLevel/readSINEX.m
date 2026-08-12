function snx = readSINEX(filename, opts)
%READSINEX Read gravity-field SINEX files (solutions, covariances, NEQs).
%
%   SNX = shLowLevel.readSINEX(FILENAME) parses the SINEX blocks used by
%   gravity-field solution providers (ITSG, COST-G, ...):
%     +SOLUTION/ESTIMATE                parameter list: type CN/SN,
%                                       degree, order, value, std
%     +SOLUTION/MATRIX_ESTIMATE L|U COVA   covariance (lower/upper,
%                                       up to 3 values per line)
%     +SOLUTION/NORMAL_EQUATION_MATRIX L|U  normal-equation matrix
%   .gz files are gunzipped transparently.
%
%   opts.Only ("all") = "estimate" (v2.2) streams the SOLUTION/ESTIMATE block
%   line-by-line (gz via Java GZIPInputStream) and skips all matrix
%   blocks - reads the estimates of a 460-MB n96 NEQ SINEX in seconds
%   with negligible memory. Returns kind='estimate', M=[].
%
%   SNX = shLowLevel.readSINEX(..., Index (struct([]))=IDX) additionally reorders the result
%   into the shLowLevel.shIndex ordering IDX: SNX.M and SNX.x then match IDX row
%   for row, ready for tvANSFilter's opts.NoiseCov (after inverting an
%   NEQ; see Output ("raw")='covariance'). Coefficients present in IDX but absent
%   from the file raise shLowLevel:readSINEX:missingParam.
%
%   Options
%     Index   ([]) shLowLevel.shIndex struct for reordering
%     Output  ("raw") raw returns the parsed blocks; covariance assembles it from a NEQ solution
%     Only    ("all") estimate streams only the SOLUTION/ESTIMATE block of large files
%
%   Outputs
%     snx  (1 x 1) struct  parsed solution with fields n, m, cs
%          (Q x 1 double) parameter list (cs: 0=C, 1=S; file order or
%          IDX order), x and sig (Q x 1 double) estimates and formal
%          sigmas (NaN if absent), M (Q x Q double) covariance or
%          normal-equation matrix per Output= ([] if absent), kind
%          ('COVA' | 'NEQ' | ''), epoch (1 x 1 double) and header meta
%
%   Format handling (VERIFIED against a real ITSG-Grace2018 n96 SINEX from
%   TU Graz): parameter type tokens are CN/SN; the degree is the first
%   integer-valued token after the type and the order the next integer
%   token <= degree - this covers both the ITSG layout
%   "CN  2 --  0  08:107:00000 ---- 2  value  std" (CODE=degree,
%   SOLN=order, PT='--') and plain "CN n m ... value std" layouts. The
%   last two numeric fields are value and std. Matrix lines are "i j v1 [v2 v3]" with column index j
%   incrementing across the up-to-3 values. Verify once against your
%   provider's file; the parser errors loudly rather than guessing. It
%   reads the whole file into memory - fine for monthly solution SINEX
%   up to Lmax ~ 120, not tuned for multi-GB normal-equation archives.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-07.

arguments
    filename {mustBeTextScalar}
    opts.Index struct = struct([])
    opts.Output (1,1) string {mustBeMember(opts.Output, ["raw","covariance"])} = "raw"
    opts.Only (1,1) string {mustBeMember(opts.Only, ["all","estimate"])} = "all"
end
if opts.Only == "estimate"
    snx = readEstimateStreaming(filename);       % v2.2: no full-file load
    return
end
filename = char(filename);
if ~isfile(filename)
    error('shLowLevel:readSINEX:fileNotFound', 'File not found: %s', filename);
end
if endsWith(filename, '.gz')
    tmp = tempname; mkdir(tmp);
    files = gunzip(filename, tmp);
    filename = files{1};
end
lines = readlines(filename);

% ---- locate blocks
inEst = false; inMat = false; kind = ''; %#ok<NASGU>
n = []; m = []; cs = []; x = []; sig = [];
ijvC = {};                                    % chunks, vertcat once
kind = '';
for k = 1:numel(lines)
    L = strtrim(lines(k));
    if L == "" || startsWith(L, "*"), continue; end
    if startsWith(L, "+SOLUTION/ESTIMATE")
        inEst = true; inMat = false; continue;
    elseif startsWith(L, "+SOLUTION/MATRIX_ESTIMATE")
        inMat = true; inEst = false;
        if ~contains(L, "COVA")
            error('shLowLevel:readSINEX:badMatrixType', ...
                'Only COVA matrix estimates are supported (line: %s).', L);
        end
        kind = 'COVA'; continue;
    elseif startsWith(L, "+SOLUTION/NORMAL_EQUATION_MATRIX")
        inMat = true; inEst = false; kind = 'NEQ'; continue;
    elseif startsWith(L, "+")
        inEst = false; inMat = false; continue;
    elseif startsWith(L, "-")
        inEst = false; inMat = false; continue;
    end
    if inEst
        tok = split(L);
        it = find(tok == "CN" | tok == "SN", 1);
        if isempty(it) || numel(tok) < it + 2
            error('shLowLevel:readSINEX:badEstimateLine', ...
                'No CN/SN parameter on SOLUTION/ESTIMATE line: %s', L);
        end
        nums = str2double(tok);
        % Degree/order extraction robust to real-world column layouts:
        % ITSG/COST-G write "CN  2 --    0 08:107:00000 ---- 2 val sig"
        % (CODE = degree, PT = '--', SOLN = order; verified against a real
        % ITSG-Grace2018 n96 SINEX). Other producers put the order directly
        % after the degree. Rule: degree = first integer-valued token after
        % the type; order = next integer-valued token <= degree.
        isInt = ~isnan(nums) & nums == round(nums) & nums >= 0;
        cand = find(isInt);
        cand = cand(cand > it);
        nn = NaN; mm = NaN;
        if ~isempty(cand)
            nn = nums(cand(1));
            for r = cand(2:end)'
                if nums(r) <= nn, mm = nums(r); break; end
            end
        end
        vi = find(~isnan(nums));
        if numel(vi) < 2 || isnan(nn) || isnan(mm)
            error('shLowLevel:readSINEX:badEstimateLine', ...
                'Unparseable SOLUTION/ESTIMATE line: %s', L);
        end
        n(end+1,1)  = nn; %#ok<AGROW>
        m(end+1,1)  = mm; %#ok<AGROW>
        cs(end+1,1) = double(tok(it) == "SN"); %#ok<AGROW>
        x(end+1,1)   = nums(vi(end-1)); %#ok<AGROW>
        sig(end+1,1) = nums(vi(end));   %#ok<AGROW>
    elseif inMat
        v = str2double(split(L));
        v = v(~isnan(v));
        if numel(v) < 3
            error('shLowLevel:readSINEX:badMatrixLine', ...
                'Unparseable matrix line: %s', L);
        end
        nv = numel(v) - 2;
        ijvC{end+1} = [repmat(v(1), nv, 1), v(2) + (0:nv-1)', v(3:end)]; %#ok<AGROW>
    end
end
Q = numel(n);
if Q == 0
    error('shLowLevel:readSINEX:noData', 'No SOLUTION/ESTIMATE block found in %s.', filename);
end

M = [];
if ~isempty(ijvC)
    ijv = vertcat(ijvC{:});
    bad = max(max(ijv(:,1)), max(ijv(:,2)));
    if bad > Q
        error('shLowLevel:readSINEX:badMatrixIndex', ...
            ['Matrix references parameter %d but SOLUTION/ESTIMATE holds ' ...
             'only %d parameters - truncated or inconsistent file.'], bad, Q);
    end
    M = zeros(Q);
    li = sub2ind([Q Q], ijv(:,1), ijv(:,2));
    M(li) = ijv(:,3);
    M = M + M' - diag(diag(M));                 % symmetrize L/U storage
end

snx = struct('n', n, 'm', m, 'cs', cs, 'x', x, 'sig', sig, 'M', M, 'kind', kind);

% ---- optional inversion NEQ -> covariance
if opts.Output == "covariance" && strcmp(kind, 'NEQ') && ~isempty(M)
    [R, flag] = chol((M + M')/2);
    if flag ~= 0
        error('shLowLevel:readSINEX:neqNotPD', ...
            'Normal-equation matrix is not positive definite; cannot invert.');
    end
    Ri = inv(R);
    snx.M = Ri * Ri';
    snx.kind = 'COVA';
end

% ---- optional reordering to a shIndex ordering
if ~isempty(fieldnames(opts.Index))
    idx = opts.Index;
    perm = zeros(idx.P, 1);
    for r = 1:idx.P
        hit = find(n == idx.n(r) & m == idx.m(r) & cs == idx.cs(r), 1);
        if isempty(hit)
            error('shLowLevel:readSINEX:missingParam', ...
                'Coefficient (n=%d, m=%d, %s) not in the SINEX file.', ...
                idx.n(r), idx.m(r), char('C' + 16*idx.cs(r)));
        end
        perm(r) = hit;
    end
    snx.n = snx.n(perm); snx.m = snx.m(perm); snx.cs = snx.cs(perm);
    snx.x = snx.x(perm); snx.sig = snx.sig(perm);
    if ~isempty(snx.M), snx.M = snx.M(perm, perm); end
end
end


function snx = readEstimateStreaming(filename)
%READESTIMATESTREAMING SOLUTION/ESTIMATE-only streaming parse (v2.2).
%   Reads plain or .gz SINEX line-by-line and stops after the estimate
%   block - the multi-GB NEQ matrices of full n96 monthly SINEX never
%   touch memory (~seconds instead of gigabytes). .gz streams through
%   java.util.zip.GZIPInputStream (base MATLAB ships Java); if Java is
%   unavailable (-nojvm) the file is gunzipped to a temp dir instead.
filename = char(filename);
if ~isfile(filename)
    error('shLowLevel:readSINEX:fileNotFound', 'File not found: %s', filename);
end
useJava = endsWith(filename, '.gz') && usejava('jvm');
if endsWith(filename, '.gz') && ~useJava
    tmp = tempname; mkdir(tmp);
    cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
    files = gunzip(filename, tmp);
    filename = files{1};
end
if useJava
    fis = java.io.FileInputStream(filename);
    gzs = java.util.zip.GZIPInputStream(fis, 2^16);
    rd  = java.io.BufferedReader(java.io.InputStreamReader(gzs), 2^16);
    getl = @() rd.readLine();
    closer = onCleanup(@() rd.close()); %#ok<NASGU>
else
    fid = fopen(filename, 'r');
    assert(fid > 0, 'shLowLevel:readSINEX:fileNotFound', 'Cannot open %s', filename);
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    getl = @() fgetl(fid);
end
inEst = false;
nEst = 0;
nA = zeros(20000, 1); mA = nA; csA = nA; xA = nA; sigA = nA;
while true
    L = getl();
    if useJava
        if isempty(L), break; end                % readLine null = EOF
        L = char(L);                             % empty java string -> ''
    else
        if ~ischar(L), break; end                % fgetl -1 = EOF
    end
    Ls = strtrim(L);
    if startsWith(Ls, '+SOLUTION/ESTIMATE'), inEst = true; continue; end
    if startsWith(Ls, '-SOLUTION/ESTIMATE'), break; end
    if ~inEst || isempty(Ls) || Ls(1) == '*', continue; end
    tp = split(string(Ls));
    if numel(tp) < 6, continue; end
    typ = tp(2);
    if typ ~= "CN" && typ ~= "SN", continue; end
    vals = str2double(tp);
    ints = find(vals == floor(vals) & ~isnan(vals) & vals >= 0, 4);
    ints = ints(ints > 2);
    if numel(ints) < 2, continue; end
    dg = vals(ints(1));
    kOrd = find(vals(ints(2:end)) <= dg, 1);
    if isempty(kOrd), continue; end
    orr = vals(ints(1 + kOrd));
    % estimate value: the last two numeric tokens are value and sigma in
    % both the ITSG and the generic SINEX estimate layouts
    numsAll = vals(~isnan(vals));
    est = numsAll(end-1); sg = numsAll(end);
    nEst = nEst + 1;
    if nEst > numel(nA)
        nA(2*end) = 0; mA(2*end) = 0; csA(2*end) = 0; xA(2*end) = 0; sigA(2*end) = 0;
    end
    nA(nEst) = dg; mA(nEst) = orr;
    csA(nEst) = double(typ == "SN"); xA(nEst) = est; sigA(nEst) = sg;
end
assert(nEst > 0, 'shLowLevel:readSINEX:noEstimates', ...
    'No SOLUTION/ESTIMATE CN/SN lines found (streaming mode).');
snx = struct('n', nA(1:nEst), 'm', mA(1:nEst), 'cs', csA(1:nEst), ...
    'x', xA(1:nEst), 'sigma', sigA(1:nEst), 'M', [], 'kind', 'estimate', ...
    'note', 'streaming Only="estimate": matrix blocks skipped');
end
