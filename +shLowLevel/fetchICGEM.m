function [file, info] = fetchICGEM(model, opts)
%FETCHICGEM Download a static gravity field model from ICGEM.
%
%   FILE = shLowLevel.fetchICGEM("EGM2008") resolves the model name against
%   shLowLevel.listICGEM (case-insensitive exact match, else unique prefix)
%   and downloads the .gfc into <dataFolder>/icgem. Existing files are
%   skipped. Load with shCoefficients.read(FILE).
%
%   FILE = shLowLevel.fetchICGEM(row) with a row of the listICGEM table skips
%   the listing round-trip.
%
%   Options
%     Dest (fullfile(shLowLevel.dataFolder(), "icgem")), Timeout (300),
%     List (table()) ([])   pass a pre-fetched listICGEM table (avoids re-listing
%                 in loops / enables the offline fixture in tests)
%     Proxy ("")  per-call proxy URL, e.g. "http://proxy:8080" (empty: MATLAB Web Preferences)
%     Quiet (false)       suppress progress output ([k/K] counter,
%                         per-file timing, failure summary)
%     Type ("static")     "static" | "temporal": which catalogue
%                         numeric selections and names resolve against
%                         (mirrors shLowLevel.listICGEM)
%     Files ("*.gfc*")    filename filter in series mode
%     Mode ("auto")       series download strategy: "auto" fetches the
%                         server's whole-series ZIP in ONE request (the
%                         rate limiter punishes hundreds of per-file
%                         requests) and falls back to file-by-file on
%                         failure; "archive" | "files" force either.
%                         Per-file is resumable with every file
%                         verified before the swap.
%     FileList (table())  series mode: use this file table (from
%                         listICGEM(Series=...)) instead of querying -
%                         for subsetting or offline mirrors
%     Pause (3)           seconds between models in bulk mode - the
%                         ICGEM server RATE-LIMITS rapid requests
%                         (HTTP 429), which looks like a stall
%     Retries (3)         retry attempts on 429/timeout, with 30/60/120
%                         second backoff
%     Update (false)  refresh existing files (safe swap: verified before replacing)
%   TIME SERIES (temporal catalogue): a catalogue row carrying the
%   'zip' column fetches a whole monthly series into
%   <dataFolder>/icgem/series/<group>_<center>_<series>/ - ready for
%   shSeries.fromFolder or shLowLevel.standardChain:
%       T  = shLowLevel.listICGEM(Type = "temporal");
%       fs = shLowLevel.fetchICGEM(T(12, :));
%       ts = shSeries.fromFolder(fileparts(fs(1)), Pattern = "*.gfc*");
%
%   Outputs
%     file  string   local path
%     info  struct: url, skipped (true if already present)
%
%   Cite the model's reference (listICGEM carries it in 'data'/ICGEM
%   pages) when publishing.
%
%   Claude (Fable 5), 2026-08-07 (v2.4.1).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    model
    opts.Dest (1,1) string = ""
    opts.Timeout (1,1) double = 300
    opts.List table = table()
    opts.Proxy (1,1) string = ""
    opts.Update (1,1) logical = false
    opts.Quiet (1,1) logical = false
    opts.Pause (1,1) double {mustBeNonnegative} = 3
    opts.Retries (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    opts.Type (1,1) string {mustBeMember(opts.Type, ...
        ["static", "temporal"])} = "static"
    opts.Files (1,1) string = "*.gfc*"
    opts.Mode (1,1) string {mustBeMember(opts.Mode, ...
        ["auto", "archive", "files"])} = "auto"
    opts.FileList table = table()
end
% ---- v3.0.0: numeric idx vector (rows of shLowLevel.listICGEM), "all", or a
% list of names fetch multiple models in one call
if isnumeric(model) || ...
        (isstring(model) && (numel(model) > 1 || model == "all"))
    if isempty(opts.List), T = shLowLevel.listICGEM(Type = opts.Type); else, T = opts.List; end
    if isnumeric(model)
        sel = model(:)';
        assert(all(sel >= 1 & sel <= height(T) & sel == round(sel)), ...
            'shLowLevel:fetchICGEM:badIdx', ...
            'Numeric selection must be catalogue rows 1..%d (see shLowLevel.listICGEM).', ...
            height(T));
        rows = T(sel, :);
    elseif isscalar(model) && model == "all"
        rows = T;
    else
        rows = T(1, []); %#ok<NASGU>
        idxs = zeros(1, numel(model));
        for k = 1:numel(model)
            h = find(strcmpi(T.name, model(k)), 1);
            assert(~isempty(h), 'shLowLevel:fetchICGEM:notFound', ...
                'Model "%s" not found (see shLowLevel.listICGEM).', model(k));
            idxs(k) = h;
        end
        rows = T(idxs, :);
    end
    K = height(rows);
    if ~opts.Quiet
        fprintf(['fetching %d ICGEM models (large fields exceed 100 MB ' ...
            'each;\npresent files are skipped, so an interrupted run ' ...
            'simply resumes)\n'], K);
    end
    fileC = repmat({strings(1, 0)}, 1, K);
    infos = cell(1, K); failed = strings(1, 0);
    for k = 1:K
        if ~opts.Quiet, fprintf('[%d/%d] ', k, K); end
        try
            [fk, infos{k}] = shLowLevel.fetchICGEM(rows(k, :), ...
                Dest = opts.Dest, Timeout = opts.Timeout, ...
                Proxy = opts.Proxy, Update = opts.Update, ...
                Quiet = opts.Quiet, Retries = opts.Retries, ...
                Type = opts.Type, Files = opts.Files, ...
                Mode = opts.Mode, List = T);
            fileC{k} = reshape(fk, 1, []);      % series return many files
            if k < K && ~infos{k}.skipped
                pause(opts.Pause);              % ICGEM rate limit
            end
        catch err
            if ismember('name', rows.Properties.VariableNames)
                nm = rows.name(k);
            elseif ismember('series', rows.Properties.VariableNames)
                nm = rows.series(k);
            else
                nm = rows.url(k);
            end
            failed(end+1) = nm + ": " + err.message; %#ok<AGROW>
            infos{k} = struct('url', rows.url(k), 'skipped', false, ...
                'updated', false, 'failed', true, 'mode', "");
            if ~opts.Quiet
                fprintf('  FAILED %s (%s)\n', nm, err.message);
            end
        end
    end
    info = [infos{:}];
    file = [fileC{:}];
    if ~opts.Quiet
        fprintf('done: %d of %d selections ok (%d files), %d failed.\n', ...
            K - numel(failed), K, numel(file), numel(failed));
        if ~isempty(failed), fprintf('  %s\n', failed); end
    end
    return
end
if istable(model)
    assert(height(model) == 1 && ismember('url', model.Properties.VariableNames), ...
        'shLowLevel:fetchICGEM:badRow', 'Table input must be ONE listICGEM row.');
    row = model;
    if ismember('zip', row.Properties.VariableNames)
        [file, info] = fetchSeries(row, opts);  % temporal catalogue row
        return
    end
else
    name = string(model);
    if isempty(opts.List)
        T = shLowLevel.listICGEM(Type = opts.Type);
    else
        T = opts.List;
    end
    if ismember('name', T.Properties.VariableNames)
        hit = find(strcmpi(T.name, name));
    else
        hit = find(strcmpi(T.series, name));
        if isempty(hit)
            hit = find(contains(T.series, name, 'IgnoreCase', true));
            if numel(hit) > 1
                error('shLowLevel:fetchICGEM:ambiguous', ...
                    '"%s" matches %d series - be more specific.', ...
                    name, numel(hit));
            end
        end
    end
    if isempty(hit)
        hit = find(startsWith(lower(T.name), lower(name)));
    end
    assert(~isempty(hit), 'shLowLevel:fetchICGEM:notFound', ...
        'Model "%s" not found at ICGEM (see shLowLevel.listICGEM).', name);
    assert(isscalar(hit), 'shLowLevel:fetchICGEM:ambiguous', ...
        '"%s" matches %d models: %s', name, numel(hit), ...
        strjoin(T.name(hit(1:min(5, end))), ', '));
    row = T(hit, :);
end
dest = opts.Dest;
if strlength(dest) == 0
    dest = string(fullfile(shLowLevel.dataFolder(), 'icgem'));
end
if ~isfolder(dest), mkdir(dest); end
[~, base, ext] = fileparts(char(row.url));
file = string(fullfile(dest, [base, ext]));
present = isfile(file);
info = struct('url', row.url, 'skipped', present && ~opts.Update, ...
    'updated', false, 'failed', false, 'mode', "model");
if present && ~opts.Update
    return
end
if ~opts.Quiet
    fprintf('  %s %s from ICGEM (large models can take minutes)...\n', ...
        ternary(present, 'updating', 'fetching'), string([base, ext]));
end
tDl = tic;
tmpf = file + ".part";
try
    fetchWithBackoff(row.url, tmpf, opts);
    if endsWith(lower(file), [".gfc", ".gfc.gz"])
        shLowLevel.shReadGFC(tmpf);                    % verify BEFORE swap
    else
        d = dir(tmpf);
        assert(~isempty(d) && d.bytes > 0, 'empty download');
    end
catch err
    if isfile(tmpf), delete(tmpf); end
    if present
        % failed refresh: the existing file stays authoritative
        warning('shLowLevel:fetchICGEM:updateFailed', '%s', err.message);
        info.skipped = true;
        return
    end
    rethrow(err);
end
movefile(tmpf, file, 'f');
info.updated = present;
info.failed = false;
if ~opts.Quiet
    d = dir(char(file));
    fprintf('  done: %.1f MB in %.0f s\n', d.bytes / 1e6, toc(tDl));
end
end

function s = ternary(tf, a, b)
if tf, s = a; else, s = b; end
end

function fetchWithBackoff(url, tmpf, opts)
% download with 429/timeout backoff; local paths copy (mirror/testing)
if isfile(url)
    copyfile(char(url), char(tmpf));
    return
end
backoff = [30, 60, 120];
attempt = 0;
while true
    try
        webFetch(url, tmpf, opts.Timeout, opts.Proxy);
        break
    catch errDl
        attempt = attempt + 1;
        retryable = contains(errDl.message, '429') || ...
            contains(errDl.message, 'Too Many Requests') || ...
            contains(lower(errDl.message), 'timed out');
        if ~retryable || attempt > opts.Retries
            rethrow(errDl);
        end
        wsec = backoff(min(attempt, numel(backoff)));
        if ~opts.Quiet
            fprintf(['  rate-limited by ICGEM - waiting %d s, ' ...
                'then retry %d/%d\n'], wsec, attempt, opts.Retries);
        end
        pause(wsec);
    end
end
end

function [files, info] = fetchSeries(row, opts)
% one temporal-catalogue row: fetch a whole monthly series into its own
% folder. Mode = "auto" (default) downloads the server's whole-series
% ZIP in ONE request - the ICGEM rate limiter punishes hundreds of
% sequential per-file requests ("too many connections" stalls observed
% in the field) - and falls back to file-by-file (resumable, each file
% verified before swap) when the archive path fails. "archive"/"files"
% force either.
dest = opts.Dest;
if strlength(dest) == 0
    tag = regexprep(char(row.group + "_" + row.center + "_" + row.series), ...
        '[^\w\-.]', '_');
    dest = string(fullfile(shLowLevel.dataFolder(), 'icgem', 'series', tag));
end
if ~isfolder(dest), mkdir(dest); end
info = struct('url', row.url, 'skipped', false, 'updated', false, ...
    'failed', false, 'mode', "");
listF = @() string(reshape({dir(fullfile(char(dest), '**', ...
    char(opts.Files))).folder}, [], 1)) + filesep + ...
    string(reshape({dir(fullfile(char(dest), '**', char(opts.Files))).name}, [], 1));
have = listF();
if ~isempty(have) && ~opts.Update
    files = reshape(have, 1, []); info.skipped = true;
    info.mode = "present";
    if ~opts.Quiet
        fprintf('  %s: %d files present, skipped\n', row.series, numel(have));
    end
    return
end
% ---- archive first: ONE request per series
if opts.Mode ~= "files"
    try
        if ~opts.Quiet
            fprintf(['  fetching series archive %s in one request ' ...
                '(can be hundreds of MB)...\n'], row.series);
        end
        tDl = tic;
        tmpf = fullfile(char(dest), 'series.zip.part');
        fetchWithBackoff(row.zip, tmpf, opts);
        unzip(tmpf, char(dest));
        delete(tmpf);
        files = reshape(listF(), 1, []);
        assert(~isempty(files), 'shLowLevel:fetchICGEM:emptySeries', ...
            'Archive of "%s" contained no files matching %s.', ...
            row.series, opts.Files);
        shLowLevel.shReadGFC(char(files(1)));   % sanity: archive intact
        info.mode = "archive";
        if ~opts.Quiet
            fprintf('  done: %d files in %.0f s\n', numel(files), toc(tDl));
        end
        return
    catch errA
        if isfile(tmpf), delete(tmpf), end
        if opts.Mode == "archive"
            rethrow(errA);
        end
        if ~opts.Quiet
            fprintf(['  archive failed (%s) - falling back to ' ...
                'file-by-file\n'], errA.message);
        end
    end
end
% ---- per-file fallback: resumable, each file verified before swap
info.mode = "files";
if isempty(opts.FileList)
    F = shLowLevel.listICGEM(Type = "temporal", Series = row.path, ...
        Timeout = opts.Timeout);
else
    F = opts.FileList;
end
pat = regexptranslate('wildcard', char(opts.Files));
keep = ~cellfun('isempty', regexp(cellstr(F.name), ['^' pat '$'], 'once'));
F = F(keep, :);
assert(height(F) > 0, 'shLowLevel:fetchICGEM:emptySeries', ...
    'Series "%s" has no files matching %s.', row.series, opts.Files);
files = strings(1, 0); nFetched = 0; nSkip = 0;
for j = 1:height(F)
    fp = fullfile(char(dest), char(F.name(j)));
    if isfile(fp) && ~opts.Update
        files(1, end+1) = string(fp); %#ok<AGROW>
        nSkip = nSkip + 1;
        continue
    end
    if ~opts.Quiet
        fprintf('  [%d/%d] %s\n', j, height(F), F.name(j));
    end
    tmpf = [fp '.part'];
    fetchedThis = true;
    try
        fetchWithBackoff(F.url(j), tmpf, opts);
        shLowLevel.shReadGFC(tmpf);             % verify BEFORE swap
        movefile(tmpf, fp, 'f');
        files(1, end+1) = string(fp); %#ok<AGROW>
        nFetched = nFetched + 1;
    catch err
        if isfile(tmpf), delete(tmpf), end
        info.failed = true;
        if ~opts.Quiet
            fprintf('  FAILED %s (%s)\n', F.name(j), err.message);
        end
        fetchedThis = false;
    end
    if j < height(F) && fetchedThis
        pause(opts.Pause);                      % ICGEM rate limit
    end
end
info.skipped = nFetched == 0 && ~info.failed;
if ~opts.Quiet
    fprintf('  series %s: %d fetched, %d present, failed: %d\n', ...
        row.series, nFetched, nSkip, info.failed);
end
end
