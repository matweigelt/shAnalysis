function model = shReadGFC(filename)
%SHREADGFC Read an ICGEM ASCII .gfc gravity field model file.
%
%   MODEL = SHREADGFC(FILENAME) parses header key:value pairs and the
%   gfc/gfct data block of a standard ICGEM-format file (GRACE/GRACE-FO
%   monthly solutions, static models, GOCE, EGM2008, etc.). Gzip-
%   compressed files (.gfc.gz) are transparently decompressed.
%
%   MODEL fields:
%     .GM, .R          gravitational parameter [m^3/s^2], reference radius [m]
%     .nmax            maximum degree found in the file
%     .C, .S           (nmax+1)x(nmax+1) coefficient matrices, C(n+1,m+1)
%                       (for gfct files: the constant/reference part only)
%     .sigmaC, .sigmaS same layout, formal errors (NaN if not given)
%     .header          struct of all raw header key:value pairs
%     .variableTerms   struct array of time-variable terms from gfct-style
%                       files (trnd/acos/asin lines), fields:
%                       .type ('trnd'/'acos'/'asin'), .n, .m, .C, .S,
%                       .t0, .t1 (validity/reference epochs, as given in
%                       the file, format/units not reinterpreted here --
%                       check the file header for the epoch convention).
%                       Empty struct array if the file is a plain (static)
%                       gfc model.
%
%   See SHEVALGFCT to evaluate C(t), S(t) at a specific epoch from a
%   gfct-type MODEL, assuming the common annual-trend + annual-periodic
%   convention (verify against your specific model's documentation --
%   ICGEM time-variable formats are not perfectly uniform across models).
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shLowLevel: Claude (Fable 5), 2026-08-07.
%   Outputs
%     model  (1,1) struct  fields: C/S/sigmaC/sigmaS (nmax+1 x nmax+1
%            double), GM/R (1,1 double), nmax (1,1 double), tide (string),
%            variable terms (trend/annual/semiannual) when present
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

readFilename = filename;
if endsWith(filename, '.gz')
    tmpDir = tempname;
    mkdir(tmpDir);
    unzipped = gunzip(filename, tmpDir);
    readFilename = unzipped{1};
    cleanupTemp = onCleanup(@() rmdir(tmpDir, 's')); %#ok<NASGU>
end

fid = fopen(readFilename, 'r');
if fid < 0
    error('shReadGFC:fileNotFound', 'Cannot open file: %s', filename);
end

header = struct();
nmax = 0;
rows = [];  % n m Cnm Snm sigC sigS  (gfc/gfct constant term)
varRows = {}; % {type, n, m, C, S, t0, t1, period}
gfctT0 = zeros(0, 3);  % ICGEM 1.0: [n m t0] reference epochs from gfct lines

% ---- v3.1.1 fast path: static files (the overwhelming majority, and
% the only ones that get LARGE - EGM2008-class bodies have millions of
% lines) are parsed in bulk. The per-line fgetl/str2double loop below
% costs ~30 s at n720 and minutes at n2190; the bulk path is ~100x
% faster and produces the identical rows matrix. Files with variable
% terms (gfct/trnd/dot/acos/asin) keep the proven line-by-line parser.
txtAll = fread(fid, inf, '*char')';
frewind(fid);
eoh = regexp(txtAll, 'end_of_head[^\n]*\n', 'end', 'once');
fastPath = ~isempty(eoh);
if fastPath                                     %#ok<ALIGN>
    % header: identical key/value semantics to the loop below
    hLines = strsplit(txtAll(1:eoh), '\n');
    for hk = 1:numel(hLines)
        tokH = strtrim(hLines{hk});
        if isempty(tokH) || startsWith(tokH, '#') || ...
                startsWith(tokH, 'end_of_head')
            continue
        end
        partsH = strsplit(tokH);
        if numel(partsH) >= 2
            keyH = matlab.lang.makeValidName(partsH{1});
            numvalH = str2double(partsH{2});
            if ~isnan(numvalH)
                header.(keyH) = numvalH;
            else
                header.(keyH) = partsH{2};
            end
        end
    end
    % body (v3.1.1): GROUP-WISE bulk parse. Lines are extracted per
    % record key (gfc/gfct/trnd/acos/asin; 'dot' is skipped exactly as
    % the line parser skips unknown keys) and each uniform group goes
    % through one sscanf - this covers variable-term files too (GRGS
    % mean fields, ~10^5..10^6 records; the per-line loop with per-line
    % calendarDuration dates took tens of minutes). Column bookkeeping:
    % V is KEY-STRIPPED, so legacy parts{k} maps to V(:, k-1). Any
    % group failing the consumed-everything check drops the whole file
    % to the legacy line parser (its skip/NaN semantics stay
    % authoritative for dirty records).
    body = txtAll(eoh+1:end);
    body = regexprep(body, '^\s*#[^\n]*\n', '', 'lineanchors');
    isV2f = isfield(header, 'format') && ...
        contains(lower(char(string(header.format))), 'icgem2.0');
    kk = {'gfc', 'gfct', 'trnd', 'acos', 'asin'};
    Gv = cell(1, 5); Gn = zeros(1, 5);
    okBulk = true;
    for ki = 1:5
        L = regexp(body, ['^[ \t]*' kk{ki} '[ \t]+[^\n]*'], ...
            'match', 'lineanchors');
        if ki == 1 && ~isempty(L)               % 'gfc\t' never hits gfct
            L = L(cellfun(@(s) isempty(regexp(s, '^[ \t]*gfct', ...
                'once')), L));
        end
        if isempty(L), continue, end
        fl = strrep(strrep(strrep(L{1}, kk{ki}, ''), 'D', 'E'), 'd', 'e');
        nc = numel(sscanf(fl, '%f'));
        J = strjoin(L, newline);
        J = strrep(J, kk{ki}, '');              % key token only
        J = strrep(strrep(J, 'D', 'E'), 'd', 'e');
        [V, ~, em, ni] = sscanf(J, '%f', [nc, Inf]);
        V = V';
        if ~isempty(em) || isempty(V) || ...
                ~all(isstrprop(J(min(ni, end):end), 'wspace'))
            okBulk = false;
            break
        end
        Gv{ki} = V; Gn(ki) = nc;
    end
    if okBulk && isempty(Gv{1}) && isempty(Gv{2})
        fclose(fid);
        error('shReadGFC:noData', ...
            'No gfc/gfct data records found in %s', filename);
    end
    varNum = zeros(0, 7);                       % n m C S t0 t1 period
    varT = cell(0, 1);
    if okBulk
        for ki = 1:5
            V = Gv{ki};
            if isempty(V) || ~okBulk, continue, end
            nc = Gn(ki); N2 = size(V, 1);
            nmax = max(nmax, max(V(:, 1)));
            if ki <= 2 && ~(ki == 2 && isV2f)
                % gfc, and gfct under ICGEM 1.0: the constant part
                r = nan(N2, 6);
                r(:, 1:min(nc, 6)) = V(:, 1:min(nc, 6));
                rows = [rows; r]; %#ok<AGROW>
                if ki == 2 && nc >= 7           % parts>=8: t0 = V(:,7)
                    gfctT0 = [gfctT0; V(:, 1:2), ...
                        shLowLevel.icgemDate2Year(V(:, 7))]; %#ok<AGROW>
                end
                continue
            end
            t0 = nan(N2, 1); t1 = t0; per = t0;
            if ki == 2                          % gfct, ICGEM 2.0 piece
                if nc < 8, okBulk = false; continue, end
                t0 = shLowLevel.icgemDate2Year(V(:, 7));
                t1 = shLowLevel.icgemDate2Year(V(:, 8));
            elseif isV2f
                % 2.0: trnd -> ... t0 t1 ; acos/asin -> ... per t0 t1
                if ki == 3
                    if nc < 8, okBulk = false; continue, end
                    t0 = shLowLevel.icgemDate2Year(V(:, 7));
                    t1 = shLowLevel.icgemDate2Year(V(:, 8));
                else
                    if nc < 9, okBulk = false; continue, end
                    per = V(:, 7);
                    t0 = shLowLevel.icgemDate2Year(V(:, 8));
                    t1 = shLowLevel.icgemDate2Year(V(:, 9));
                end
            else
                % 1.0 disambiguation (documented in the line parser):
                % parts>=9 (nc>=8): t0/t1 = V(:,7:8), period 1 yr for
                % acos/asin; parts==8 & ~trnd (nc==7): V(:,7) = period
                if nc >= 8
                    t0 = shLowLevel.icgemDate2Year(V(:, 7));
                    t1 = shLowLevel.icgemDate2Year(V(:, 8));
                    if ki > 3, per = ones(N2, 1); end
                elseif nc == 7 && ki > 3
                    per = V(:, 7);
                end
            end
            varNum = [varNum; V(:, 1:2), V(:, 3:4), t0, t1, per]; %#ok<AGROW>
            varT = [varT; repmat(kk(ki), N2, 1)]; %#ok<AGROW>
        end
    end
    if okBulk
        if ~isempty(varNum)
            varRows = [varT, num2cell(varNum)];
        end
        fclose(fid);
    else
        fastPath = false;                       % dirty -> line parser
        header = struct(); rows = []; nmax = 0;
        varRows = {}; gfctT0 = zeros(0, 3);
        varNum = zeros(0, 7); varT = cell(0, 1); %#ok<NASGU>
    end
end
if ~fastPath
line = fgetl(fid);
inHeader = true;
while ischar(line)
    tok = strtrim(line);
    if inHeader
        if startsWith(tok, 'end_of_head')
            inHeader = false;
        elseif ~isempty(tok) && ~startsWith(tok, '#')
            parts = strsplit(tok);
            if numel(parts) >= 2
                key = matlab.lang.makeValidName(parts{1});
                val = parts{2};
                numval = str2double(val);
                if ~isnan(numval)
                    header.(key) = numval;
                else
                    header.(key) = val;
                end
            end
        end
    else
        if ~isempty(tok) && ~startsWith(tok, '#')
            parts = strsplit(tok);
            key = lower(parts{1});
            % FORTRAN D-exponents: str2double returns NaN for them (CI-
            % verified), which silently corrupted D-files to NaN above
            % the switch degree. Normalize numeric parts only (the key
            % 'dot' must survive).
            parts(2:end) = strrep(strrep(parts(2:end), 'D', 'E'), 'd', 'e');
            isV2 = isfield(header, 'format') && ischar(header.format) ...
                && contains(lower(header.format), 'icgem2.0');
            if strcmp(key, 'gfc') || (strcmp(key, 'gfct') && ~isV2)
                n = str2double(parts{2});
                m = str2double(parts{3});
                Cnm = str2double(parts{4});
                Snm = str2double(parts{5});
                if numel(parts) >= 7
                    sC = str2double(parts{6});
                    sS = str2double(parts{7});
                else
                    sC = NaN; sS = NaN;
                end
                rows(end+1, :) = [n m Cnm Snm sC sS]; %#ok<AGROW>
                nmax = max(nmax, n);
                if strcmp(key, 'gfct') && numel(parts) >= 8
                    % ICGEM 1.0: reference epoch of this coefficient
                    gfctT0(end+1, :) = ...
                        [n m shLowLevel.icgemDate2Year(str2double(parts{8}))]; %#ok<AGROW>
                end
            elseif strcmp(key, 'gfct')                     % ICGEM 2.0 piece
                n = str2double(parts{2});
                m = str2double(parts{3});
                t0 = shLowLevel.icgemDate2Year(str2double(parts{8}));
                t1 = shLowLevel.icgemDate2Year(str2double(parts{9}));
                varRows(end+1,:) = {'gfct', n, m, str2double(parts{4}), ...
                    str2double(parts{5}), t0, t1, NaN}; %#ok<AGROW>
                nmax = max(nmax, n);
            elseif strcmp(key, 'trnd') || strcmp(key, 'acos') || strcmp(key, 'asin')
                n = str2double(parts{2});
                m = str2double(parts{3});
                Cval = str2double(parts{4});
                Sval = str2double(parts{5});
                t0 = NaN; t1 = NaN; period = NaN;
                if isV2
                    % 2.0: trnd -> ... sigC sigS t0 t1
                    %      acos/asin -> ... sigC sigS period[yr] t0 t1
                    if strcmp(key, 'trnd')
                        t0 = shLowLevel.icgemDate2Year(str2double(parts{8}));
                        t1 = shLowLevel.icgemDate2Year(str2double(parts{9}));
                    else
                        period = str2double(parts{8});
                        t0 = shLowLevel.icgemDate2Year(str2double(parts{9}));
                        t1 = shLowLevel.icgemDate2Year(str2double(parts{10}));
                    end
                else
                    % 1.0 conventions vary; disambiguation (documented):
                    %  - 8 columns on acos/asin: col 8 is the PERIOD [yr],
                    %    reference epoch comes from the gfct line (EIGEN
                    %    style)
                    %  - >= 9 columns: cols 8/9 are t0/t1 (period = 1 yr),
                    %    preserving pre-v2.1 behavior for existing files
                    if numel(parts) >= 9
                        t0 = shLowLevel.icgemDate2Year(str2double(parts{8}));
                        t1 = shLowLevel.icgemDate2Year(str2double(parts{9}));
                        if ~strcmp(key, 'trnd'), period = 1.0; end
                    elseif numel(parts) == 8 && ~strcmp(key, 'trnd')
                        period = str2double(parts{8});
                    end
                end
                varRows(end+1,:) = {key, n, m, Cval, Sval, t0, t1, period}; %#ok<AGROW>
                nmax = max(nmax, n);
            end
        end
    end
    line = fgetl(fid);
end
fclose(fid);
end                                             % fastPath / legacy

if isempty(rows) && isempty(varRows)
    error('shReadGFC:noData', 'No gfc/gfct data records found in %s', filename);
end

C = zeros(nmax+1, nmax+1);
S = zeros(nmax+1, nmax+1);
sigmaC = nan(nmax+1, nmax+1);
sigmaS = nan(nmax+1, nmax+1);
for k = 1:size(rows,1)
    n = rows(k,1); m = rows(k,2);
    C(n+1,m+1) = rows(k,3);
    S(n+1,m+1) = rows(k,4);
    sigmaC(n+1,m+1) = rows(k,5);
    sigmaS(n+1,m+1) = rows(k,6);
end

model.header = header;
model.nmax = nmax;
model.C = C;
model.S = S;
model.sigmaC = sigmaC;
model.sigmaS = sigmaS;

if isempty(varRows)
    model.variableTerms = struct('type',{},'n',{},'m',{},'C',{},'S',{}, ...
        't0',{},'t1',{},'period',{});
else
    % vectorized assembly (v3.1.1): the per-element struct growth was
    % O(N^2) at GRGS scale. Patch missing t0 from the 1.0 gfct epochs
    % first, then build the struct array in one call.
    vN = cell2mat(varRows(:, 2:8));
    if ~isempty(gfctT0)
        miss = isnan(vN(:, 5));
        if any(miss)
            [tf, loc] = ismember(vN(miss, 1:2), gfctT0(:, 1:2), 'rows');
            idxm = find(miss);
            vN(idxm(tf), 5) = gfctT0(loc(tf), 3);
        end
    end
    vt = struct('type', varRows(:, 1)', ...
        'n', num2cell(vN(:, 1))', 'm', num2cell(vN(:, 2))', ...
        'C', num2cell(vN(:, 3))', 'S', num2cell(vN(:, 4))', ...
        't0', num2cell(vN(:, 5))', 't1', num2cell(vN(:, 6))', ...
        'period', num2cell(vN(:, 7))');
    model.variableTerms = vt;

    % ICGEM 2.0: gfct pieces live in variableTerms; expose the most recent
    % piece in C/S so static use of the model remains meaningful
    isPiece = strcmp({vt.type}, 'gfct');
    if any(isPiece)
        pcs = vt(isPiece);
        [~, ord] = sort([pcs.t0]);
        for q = ord(:)'
            C(pcs(q).n+1, pcs(q).m+1) = pcs(q).C;
            S(pcs(q).n+1, pcs(q).m+1) = pcs(q).S;
        end
        model.C = C; model.S = S;
    end
end

if isfield(header, 'earth_gravity_constant')
    model.GM = header.earth_gravity_constant;
elseif isfield(header, 'GM')
    model.GM = header.GM;
else
    model.GM = 3.986004415e14; % fallback, GRACE-era convention
end

if isfield(header, 'radius')
    model.R = header.radius;
elseif isfield(header, 'R')
    model.R = header.R;
else
    model.R = 6378136.3; % fallback, GRACE-era convention
end

if isfield(header, 'x_modelname') || isfield(header, 'modelname')
    if isfield(header,'modelname'), model.modelname = header.modelname;
    else, model.modelname = header.x_modelname; end
else
    model.modelname = '';
end

end
