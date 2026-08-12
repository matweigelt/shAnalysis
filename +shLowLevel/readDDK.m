function W = readDDK(source, opts)
%READDDK Load a DDK (Kusche) anisotropic decorrelation filter.
%
%   W = shLowLevel.readDDK(SOURCE) loads a DDK filter into the toolbox's
%   order-block container, ready for shLowLevel.applyDDK / g.applyDDK /
%   ts.applyDDK. The DDK filters (Kusche 2007; Kusche et al. 2009,
%   DDK1..DDK8) are the standard published anisotropic alternative to
%   Gaussian smoothing and the usual benchmark against tvANS.
%
%   Accepted SOURCE forms:
%     binary Wbd file      the released filter files (Rietbroek/Kusche
%                          'BIN' format, type BDFULLV0), e.g.
%                          Wbd_2-120.a_1d12p_4 from
%                          github.com/strawpants/GRACE-filter (MIT).
%                          Parsed natively (v2.2.1): little/big endian
%                          autodetect, versions 2.1 and >= 2.4; verified
%                          against the repository-documented pack values
%                          of Wbd_2-120.a_1d12p_4 to 4.4e-16.
%                          DDK number -> file:  DDK1=a_1d14p_4,
%                          DDK2=1d13, DDK3=1d12, DDK4=5d11, DDK5=1d11,
%                          DDK6=5d10, DDK7=1d10, DDK8=5d9.
%     .mat file / struct   container layout below (converted once via
%                          W = shLowLevel.readDDK('Wbd_...'); save ddk3.mat W)
%     ASCII file           the documented plain-text exchange layout:
%                            # header lines (ignored)
%                            block m cs n1 n2      (order, 0=C/1=S,
%                                                   degree range)
%                            <(n2-n1+1) rows of (n2-n1+1) values>
%                          repeated per block.
%
%   Options
%     Nmax ([]) truncate the filter to degrees <= Nmax (leading
%               submatrix per block - identical to the reference
%               filterSH behavior for inputs with nmax below the
%               filter's Lmax). Required to apply an Lmax-120 DDK to
%               n96 series without an nmaxMismatch error.
%
%   Container layout (returned / expected in .mat):
%     W.nmax    (1,1) double
%     W.blocks  struct array: m, cs, n (row vector of degrees), M
%               (square filter matrix over those degrees, applied as
%               cFiltered(n) = M * c(n) within the order/cs block)
%     W.name    string
%
%   Degrees below the filter's Lmin (typically 2) are not stored and pass
%   through unfiltered - consistent with the reference implementation.
%
%   Claude (Fable 5), 2026-08-07 (binary support v2.2.1, same day).
%   Inputs
%     source          DDK filter file or folder the coefficients are
%                     read from
%
%   Outputs
%     W          struct: nmax (1 x 1), name, blocks (1 x 241) struct array with m, cs, n (degrees), M (square block matrix)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    source
    opts.Nmax double = []
end
if (ischar(source) || isstring(source)) && ...
        ~isempty(regexp(char(source), '^[Dd][Dd][Kk][1-8]$', 'once'))
    % name form "DDK<n>" (v2.4.1): resolve against <dataFolder>/DDK and
    % the shipped tests/test_data (DDK3), else point at shLowLevel.fetchDDK
    sc = char(source);
    k = str2double(sc(4));
    fn = shLowLevel.ddkNames(); fn = fn(k);
    cand = [ ...
        string(fullfile(shLowLevel.dataFolder(), 'DDK', fn))
        string(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
            'tests', 'test_data', fn))];
    hit = cand(arrayfun(@isfile, cand));
    assert(~isempty(hit), 'shLowLevel:readDDK:notFetched', ...
        ['%s (%s) is not in %s or the shipped test data.\n' ...
         'Download it once with shLowLevel.fetchDDK(%d).'], ...
        upper(char(source)), fn, cand(1), k);
    W = shLowLevel.readDDK(hit(1), Nmax = opts.Nmax);
    return
end
if isstruct(source)
    W = source;
elseif (ischar(source) || isstring(source)) && endsWith(source, '.mat')
    assert(isfile(source), 'shLowLevel:readDDK:fileNotFound', 'File not found: %s', source);
    S = load(source);
    fn = fieldnames(S);
    assert(numel(fn) >= 1, 'shLowLevel:readDDK:emptyMat', 'MAT file holds no variables.');
    W = S.(fn{1});
elseif (ischar(source) || isstring(source)) && isWbdBinary(source)
    W = readWbdBinary(char(source));
elseif ischar(source) || isstring(source)
    assert(isfile(source), 'shLowLevel:readDDK:fileNotFound', 'File not found: %s', source);
    lines = readlines(source);
    W = struct('nmax', 0, 'blocks', struct('m', {}, 'cs', {}, 'n', {}, 'M', {}), ...
        'name', string(source));
    k = 1;
    while k <= numel(lines)
        L = strtrim(lines(k)); k = k + 1;
        if L == "" || startsWith(L, "#"), continue; end
        p = split(L);
        if p(1) ~= "block" || numel(p) < 5
            error('shLowLevel:readDDK:badLine', 'Expected "block m cs n1 n2", got: %s', L);
        end
        m = str2double(p(2)); cs = str2double(p(3));
        n1 = str2double(p(4)); n2 = str2double(p(5));
        nd = n2 - n1 + 1;
        M = zeros(nd);
        for r = 1:nd
            v = str2double(split(strtrim(lines(k)))); k = k + 1;
            v = v(~isnan(v));
            assert(numel(v) == nd, 'shLowLevel:readDDK:badMatrixRow', ...
                'Block (m=%d,cs=%d): row %d has %d of %d values.', m, cs, r, numel(v), nd);
            M(r, :) = v;
        end
        W.blocks(end+1) = struct('m', m, 'cs', cs, 'n', n1:n2, 'M', M); %#ok<AGROW>
        W.nmax = max(W.nmax, n2);
    end
    assert(~isempty(W.blocks), 'shLowLevel:readDDK:noData', 'No blocks found in %s.', source);
else
    error('shLowLevel:readDDK:badSource', 'SOURCE must be a struct, .mat file, or ASCII file.');
end
assert(isfield(W, 'blocks') && isfield(W, 'nmax'), 'shLowLevel:readDDK:badContainer', ...
    'W must carry fields nmax and blocks(m, cs, n, M).');
if ~isempty(opts.Nmax)
    keepB = true(1, numel(W.blocks));
    for k = 1:numel(W.blocks)
        keep = W.blocks(k).n <= opts.Nmax;
        if ~any(keep)
            keepB(k) = false;
            continue;
        end
        W.blocks(k).n = W.blocks(k).n(keep);
        W.blocks(k).M = W.blocks(k).M(keep, keep);   % leading submatrix
    end
    W.blocks = W.blocks(keepB);
    W.nmax = min(W.nmax, opts.Nmax);
end
end

function tf = isWbdBinary(filename)
%ISWBDBINARY True if the file starts with the BIN endian marker 18754.
tf = false;
if ~isfile(filename), return; end
fid = fopen(filename, 'r', 'ieee-le');
if fid < 0, return; end
marker = fread(fid, 1, 'uint16');
fclose(fid);
tf = isequal(marker, 18754) || isequal(marker, 16969); % LE or byte-swapped
end

function W = readWbdBinary(filename)
%READWBDBINARY Native parser for the Rietbroek/Kusche BIN BDFULLV0 files.
%   Layout learned from the MIT-licensed reference reader (read_BIN.m,
%   github.com/strawpants/GRACE-filter); implementation original.
%   Python-validated on the real Wbd_2-120.a_1d12p_4 against the
%   repository-documented first pack values (4.4e-16) with all 241
%   blocks structurally consistent and the pack vector fully consumed.
fid = fopen(filename, 'r', 'ieee-le');
assert(fid > 0, 'shLowLevel:readDDK:fileNotFound', 'Cannot open %s', filename);
closer = onCleanup(@() fclose(fid));
marker = fread(fid, 1, 'uint16');
if marker ~= 18754
    clear closer;                                 % reopen big-endian
    fid = fopen(filename, 'r', 'ieee-be');
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    marker = fread(fid, 1, 'uint16');
    assert(isequal(marker, 18754), 'shLowLevel:readDDK:badBinary', ...
        'Endian marker mismatch in %s.', filename);
end
version = ['BI', fread(fid, 6, 'uint8=>char')'];
ver = str2double(version(5:8));
mtype = fread(fid, 8, 'uint8=>char')';
fread(fid, 80, 'uint8');                          % descr
mi = fread(fid, 4, 'uint32');
nints = mi(1); ndbls = mi(2); nval1 = mi(3); nval2 = mi(4);
if ver < 2.4
    mi = fread(fid, 2, 'uint32');
else
    mi = fread(fid, 2, 'uint64');
end
pval1 = mi(1);
if ver <= 2.1
    nvec = 0; nread = 0;
    nval2 = nval1; %#ok<NASGU>
else
    nvec = fread(fid, 1, 'int32');
    nread = fread(fid, 1, 'int32');
end
assert(strncmp(mtype, 'BD', 2), 'shLowLevel:readDDK:badBinary', ...
    'Expected a block-diagonal (BD*) matrix, got type "%s".', mtype);
nblocks = fread(fid, 1, 'int32');
if nread > 0, fread(fid, 80 * nread, 'uint8'); end
Lmin = 2; Lmax = [];
if nints > 0
    intsD = reshape(fread(fid, 24 * nints, 'uint8=>char'), 24, nints)';
    if ver <= 2.4
        ints = fread(fid, nints, 'int32');
    else
        ints = fread(fid, nints, 'int64');
    end
    for i = 1:nints
        if strncmp(intsD(i, :), 'Lmax', 4), Lmax = ints(i); end
        if strncmp(intsD(i, :), 'Lmin', 4), Lmin = ints(i); end
    end
end
assert(~isempty(Lmax), 'shLowLevel:readDDK:badBinary', 'No Lmax metadata found.');
if ndbls > 0
    fread(fid, 24 * ndbls, 'uint8');
    fread(fid, ndbls, 'double');
end
fread(fid, 24 * nval1, 'uint8');                  % side1_d
blockind = fread(fid, nblocks, 'int32');
if any(strcmp(mtype, {'BDFULLV0', 'BDFULLVN'})) && ver > 2.2
    fread(fid, 24 * nval2, 'uint8');              % side2_d
end
for i = 1:nvec, fread(fid, nval1, 'double'); end
pack = fread(fid, pval1, 'double');
assert(numel(pack) == pval1, 'shLowLevel:readDDK:badBinary', ...
    'Truncated pack section (%d of %d doubles).', numel(pack), pval1);

W = struct('nmax', Lmax, ...
    'blocks', struct('m', {}, 'cs', {}, 'n', {}, 'M', {}), ...
    'name', string(filename));
last = 0; idx0 = 0;
for iblk = 1:nblocks
    m = floor(iblk / 2);
    if iblk == 1
        cs = 0;                                   % C, m = 0
    elseif mod(iblk, 2) == 0
        cs = 0;                                   % C
    else
        cs = 1;                                   % S
    end
    sz = blockind(iblk) - last;
    nmin = max(Lmin, m);
    assert(sz == Lmax - nmin + 1, 'shLowLevel:readDDK:badBinary', ...
        'Block %d size %d inconsistent with degrees %d..%d.', iblk, sz, nmin, Lmax);
    M = reshape(pack(idx0 + 1: idx0 + sz^2), sz, sz);  % column-major = fread order
    W.blocks(end+1) = struct('m', m, 'cs', cs, 'n', nmin:Lmax, 'M', M); %#ok<AGROW>
    last = blockind(iblk); idx0 = idx0 + sz^2;
end
assert(idx0 == pval1, 'shLowLevel:readDDK:badBinary', ...
    'Pack vector not fully consumed (%d of %d).', idx0, pval1);
end
