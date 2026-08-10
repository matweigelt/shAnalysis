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
%     Quiet (false)       suppress progress output (sizes, [k/K]
%                         counter, per-file timing, failure summary)
%     Update (false)  refresh existing files (safe swap: verified before replacing)
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
end
% ---- v3.0.0: numeric idx vector (rows of shLowLevel.listICGEM), "all", or a
% list of names fetch multiple models in one call
if isnumeric(model) || ...
        (isstring(model) && (numel(model) > 1 || model == "all"))
    if isempty(opts.List), T = shLowLevel.listICGEM(); else, T = opts.List; end
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
    file = strings(1, K); infos = cell(1, K); failed = strings(1, 0);
    for k = 1:K
        if ~opts.Quiet, fprintf('[%d/%d] ', k, K); end
        try
            [file(k), infos{k}] = shLowLevel.fetchICGEM(rows(k, :), ...
                Dest = opts.Dest, Timeout = opts.Timeout, ...
                Proxy = opts.Proxy, Update = opts.Update, ...
                Quiet = opts.Quiet, List = T);
        catch err
            if ismember('name', rows.Properties.VariableNames)
                nm = rows.name(k);
            else
                nm = rows.url(k);
            end
            failed(end+1) = nm + ": " + err.message; %#ok<AGROW>
            infos{k} = struct('url', rows.url(k), 'skipped', false, ...
                'updated', false, 'failed', true);
            if ~opts.Quiet
                fprintf('  FAILED %s (%s)\n', nm, err.message);
            end
        end
    end
    info = [infos{:}];
    keep = strlength(file) > 0;
    file = file(keep);
    if ~opts.Quiet
        fprintf('done: %d ok, %d failed.\n', nnz(keep), numel(failed));
        if ~isempty(failed), fprintf('  %s\n', failed); end
    end
    return
end
if istable(model)
    assert(height(model) == 1 && ismember('url', model.Properties.VariableNames), ...
        'shLowLevel:fetchICGEM:badRow', 'Table input must be ONE listICGEM row.');
    row = model;
else
    name = string(model);
    if isempty(opts.List)
        T = shLowLevel.listICGEM();
    else
        T = opts.List;
    end
    hit = find(strcmpi(T.name, name));
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
    'updated', false, 'failed', false);
if present && ~opts.Update
    return
end
if ~opts.Quiet
    nb = headBytes(row.url, min(opts.Timeout, 15), opts.Proxy);
    sz = '';
    if isfinite(nb), sz = sprintf(' (%.1f MB)', nb / 1e6); end
    fprintf('  %s %s%s from ICGEM...\n', ...
        ternary(present, 'updating', 'fetching'), string([base, ext]), sz);
end
tDl = tic;
tmpf = file + ".part";
try
    webFetch(row.url, tmpf, opts.Timeout, opts.Proxy);
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
