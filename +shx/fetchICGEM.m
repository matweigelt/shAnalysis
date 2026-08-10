function [file, info] = fetchICGEM(model, opts)
%FETCHICGEM Download a static gravity field model from ICGEM.
%
%   FILE = shx.fetchICGEM("EGM2008") resolves the model name against
%   shx.listICGEM (case-insensitive exact match, else unique prefix)
%   and downloads the .gfc into <dataFolder>/icgem. Existing files are
%   skipped. Load with shCoefficients.read(FILE).
%
%   FILE = shx.fetchICGEM(row) with a row of the listICGEM table skips
%   the listing round-trip.
%
%   Options
%     Dest (fullfile(shx.dataFolder(), "icgem")), Timeout (300),
%     List ([])   pass a pre-fetched listICGEM table (avoids re-listing
%                 in loops / enables the offline fixture in tests)
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
end
if istable(model)
    assert(height(model) == 1 && ismember('url', model.Properties.VariableNames), ...
        'shx:fetchICGEM:badRow', 'Table input must be ONE listICGEM row.');
    row = model;
else
    name = string(model);
    if isempty(opts.List)
        T = shx.listICGEM();
    else
        T = opts.List;
    end
    hit = find(strcmpi(T.name, name));
    if isempty(hit)
        hit = find(startsWith(lower(T.name), lower(name)));
    end
    assert(~isempty(hit), 'shx:fetchICGEM:notFound', ...
        'Model "%s" not found at ICGEM (see shx.listICGEM).', name);
    assert(isscalar(hit), 'shx:fetchICGEM:ambiguous', ...
        '"%s" matches %d models: %s', name, numel(hit), ...
        strjoin(T.name(hit(1:min(5, end))), ', '));
    row = T(hit, :);
end
dest = opts.Dest;
if strlength(dest) == 0
    dest = string(fullfile(shx.dataFolder(), 'icgem'));
end
if ~isfolder(dest), mkdir(dest); end
[~, base, ext] = fileparts(char(row.url));
file = string(fullfile(dest, [base, ext]));
present = isfile(file);
info = struct('url', row.url, 'skipped', present && ~opts.Update, ...
    'updated', false);
if present && ~opts.Update
    return
end
fprintf('  %s %s from ICGEM...\n', ...
    ternary(present, 'updating', 'fetching'), string([base, ext]));
tmpf = file + ".part";
try
    webFetch(row.url, tmpf, opts.Timeout, opts.Proxy);
    if endsWith(lower(file), [".gfc", ".gfc.gz"])
        shx.shReadGFC(tmpf);                    % verify BEFORE swap
    else
        d = dir(tmpf);
        assert(~isempty(d) && d.bytes > 0, 'empty download');
    end
catch err
    if isfile(tmpf), delete(tmpf); end
    if present
        % failed refresh: the existing file stays authoritative
        warning('shx:fetchICGEM:updateFailed', '%s', err.message);
        info.skipped = true;
        return
    end
    rethrow(err);
end
movefile(tmpf, file, 'f');
info.updated = present;
end

function s = ternary(tf, a, b)
if tf, s = a; else, s = b; end
end
