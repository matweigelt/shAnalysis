function [files, info] = fetchDDK(which, opts)
%FETCHDDK Download DDK filter matrices (Wbd binaries) on demand.
%
%   FILES = shx.fetchDDK(1:8) fetches the released anisotropic DDK
%   decorrelation filters DDK1..DDK8 (Kusche 2007; Kusche et al. 2009;
%   ~9.8 MB each) from the MIT-licensed strawpants/GRACE-filter
%   repository into <dataFolder>/DDK. Existing files are skipped unless
%   Update=true, which re-downloads them with a safe swap (the fresh
%   file is parse-verified by shx.readDDK before replacing the old one).
%   DDK3 additionally ships inside the toolbox (tests/test_data) and
%   never needs fetching.
%
%   Ordering strong -> weak smoothing, regularization a*l^4:
%     DDK1 1e14   DDK2 1e13   DDK3 1e12   DDK4 5e11
%     DDK5 1e11   DDK6 5e10   DDK7 1e10   DDK8 5e9
%
%   Inputs
%     which  (1,:) double  subset of 1..8
%   Options
%     Dest (fullfile(shx.dataFolder(), "DDK")), Timeout (120), Quiet
%   Outputs
%     files  (1,:) string  present file paths (new + existing)
%     info   struct: fetched, updated, skipped, names, url
%
%   Load with W = shx.readDDK("DDK5") (resolves against Dest) or with
%   the explicit path. Cite Kusche et al. when publishing.
%
%   Claude (Fable 5), 2026-08-07 (v2.4.1).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    which (1,:) double {mustBeInteger, mustBeInRange(which, 1, 8)}
    opts.Dest (1,1) string = ""
    opts.Timeout (1,1) double = 120
    opts.Proxy (1,1) string = ""
    opts.Update (1,1) logical = false
    opts.Quiet (1,1) logical = false
end
names = shx.ddkNames();
dest = opts.Dest;
if strlength(dest) == 0
    dest = string(fullfile(shx.dataFolder(), 'DDK'));
end
if ~isfolder(dest), mkdir(dest); end
base = "https://raw.githubusercontent.com/strawpants/GRACE-filter/master/data/DDK";
wo = weboptions('Timeout', opts.Timeout);
files = strings(1, 0); fetched = files; skipped = files; updated = files;
for k = unique(which(:)')
    fn = names(k);
    fp = fullfile(dest, fn);
    present = isfile(fp);
    if present && ~opts.Update
        skipped(end+1) = string(fp); files(end+1) = string(fp); %#ok<AGROW>
        continue
    end
    if ~opts.Quiet
        fprintf('  %s DDK%d (%s, ~10 MB)...\n', ...
            ternary(present, 'updating', 'fetching'), k, fn);
    end
    tmpf = fp + ".part";
    try
        webFetch(base + "/" + fn, tmpf, opts.Timeout, opts.Proxy);
        shx.readDDK(tmpf);                      % verify BEFORE swap
    catch err
        if isfile(tmpf), delete(tmpf); end
        if present
            % failed refresh: the existing file stays authoritative
            warning('shx:fetchDDK:updateFailed', 'DDK%d: %s', k, err.message);
            skipped(end+1) = string(fp); files(end+1) = string(fp); %#ok<AGROW>
            continue
        end
        rethrow(err);                           % new file: surface the error
    end
    movefile(tmpf, fp, 'f');
    files(end+1) = string(fp); %#ok<AGROW>
    if present
        updated(end+1) = string(fp); %#ok<AGROW>
    else
        fetched(end+1) = string(fp); %#ok<AGROW>
    end
end
info = struct('fetched', fetched, 'updated', updated, 'skipped', skipped, ...
    'names', names(unique(which(:)')), 'url', base, 'dest', dest);
end

function s = ternary(tf, a, b)
if tf, s = a; else, s = b; end
end
