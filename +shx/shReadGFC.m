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
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shx: Claude (Fable 5), 2026-08-07.
%   Outputs
%     C          (nmax+1 x nmax+1) double   cosine coefficients (static part)
%     S          (nmax+1 x nmax+1) double   sine coefficients
%     header     struct: GM, R, nmax, tide system, errors flag, variableTerms for gfct files
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
                        [n m shx.icgemDate2Year(str2double(parts{8}))]; %#ok<AGROW>
                end
            elseif strcmp(key, 'gfct')                     % ICGEM 2.0 piece
                n = str2double(parts{2});
                m = str2double(parts{3});
                t0 = shx.icgemDate2Year(str2double(parts{8}));
                t1 = shx.icgemDate2Year(str2double(parts{9}));
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
                        t0 = shx.icgemDate2Year(str2double(parts{8}));
                        t1 = shx.icgemDate2Year(str2double(parts{9}));
                    else
                        period = str2double(parts{8});
                        t0 = shx.icgemDate2Year(str2double(parts{9}));
                        t1 = shx.icgemDate2Year(str2double(parts{10}));
                    end
                else
                    % 1.0 conventions vary; disambiguation (documented):
                    %  - 8 columns on acos/asin: col 8 is the PERIOD [yr],
                    %    reference epoch comes from the gfct line (EIGEN
                    %    style)
                    %  - >= 9 columns: cols 8/9 are t0/t1 (period = 1 yr),
                    %    preserving pre-v2.1 behavior for existing files
                    if numel(parts) >= 9
                        t0 = shx.icgemDate2Year(str2double(parts{8}));
                        t1 = shx.icgemDate2Year(str2double(parts{9}));
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
    vt = struct('type',{},'n',{},'m',{},'C',{},'S',{},'t0',{},'t1',{}, ...
        'period',{});
    for k = 1:size(varRows,1)
        vt(k).type = varRows{k,1};
        vt(k).n = varRows{k,2};
        vt(k).m = varRows{k,3};
        vt(k).C = varRows{k,4};
        vt(k).S = varRows{k,5};
        vt(k).t0 = varRows{k,6};
        vt(k).t1 = varRows{k,7};
        vt(k).period = varRows{k,8};
        if isnan(vt(k).t0) && ~isempty(gfctT0)
            hit = gfctT0(:,1) == vt(k).n & gfctT0(:,2) == vt(k).m;
            if any(hit)
                vt(k).t0 = gfctT0(find(hit, 1), 3);
            end
        end
    end
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
