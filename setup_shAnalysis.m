function summary = setup_shAnalysis(opts)
%SETUP_SHANALYSIS One-call setup: path + staged data downloads.
%
%   setup_shAnalysis                        adds the toolbox to the path
%                                           for THIS session only
%   setup_shAnalysis(Permanent (false) = true)      also persists it (savepath)
%   setup_shAnalysis(Download ("none") = "core")     + TN-14 and TN-13 files
%   setup_shAnalysis(Download = "filters")  + DDK (3) filter matrices
%   setup_shAnalysis(Download = "starter")  + ITSG starter months
%   setup_shAnalysis(DryRun (false) = true, ...)    print/return the plan, do
%                                           NOTHING (no path change, no
%                                           folders, no network)
%
%   Download levels are CUMULATIVE ("several levels"):
%     "none"     path only (default).
%     "core"     shLowLevel.fetchTN: GSFC TN-14 (C20/C30) + TN-13 geocenter
%                files (GFZ/CSR/JPL, ~250 KB total) - the low-degree
%                completion chain. Every file is verified by parse.
%     "filters"  core + shLowLevel.fetchDDK(DDK): anisotropic DDK
%                decorrelation matrices (~10 MB each).
%     "starter"  filters + shLowLevel.fetchITSG(Months (["2008-04", "2025-12"]), Nmax (96)=Nmax): a small
%                monthly-solution starter set - by default the n96
%                companions of the two shipped n60 fixture months
%                (2008-04 GRACE, 2025-12 GRACE-FO).
%   Downloads land in shLowLevel.dataFolder() subfolders (TN/, DDK/,
%   itsg_series/); existing files are skipped, so re-running is cheap
%   and idempotent. Failures are collected per file and reported - an
%   offline machine still gets a working path setup plus the shipped
%   fixtures.
%
%   Options
%     Permanent (false)   addpath + savepath; if savepath fails (no
%                         write permission on pathdef.m) a warning
%                         explains the startup.m alternative
%     Download ("none")   "none" | "core" | "filters" | "starter"
%     DDK (3)             subset of 1..8 for the "filters" level
%     Providers (["GFZ","CSR","JPL"])  TN-13 providers for "core"
%     Months (["2008-04","2025-12"])   months for "starter"
%     Nmax (96)           degree of the ITSG starter files (60 or 96)
%     Docs (false)        builddocsearchdb on html/ so the toolbox
%                         pages appear in the MATLAB help-browser search
%     DryRun (false)      report the plan only - zero side effects
%     Update (false)     forward to the fetchers: refresh existing files
%                       with a safe, parse-verified swap (temporal files
%                       like TN-13/TN-14 grow monthly)
%     SeriesFolder ("")   folder of monthly .gfc files, exported as the
%                         SHX_SERIES_FOLDER environment variable AND
%                         persisted via setpref('shAnalysis',
%                         'SeriesFolder', ...) - preferences survive
%                         MATLAB restarts, no startup.m needed. This
%                         enables the opt-in trend regression in
%                         testScience: a trend needs a real monthly
%                         series, which is far too large to ship as a
%                         fixture, and the fixture suite cannot see the
%                         class of failure that only appears on a real
%                         series of real length. Deliberately NOT wired
%                         to the persistent data folder, so a routine
%                         acceptance run never touches a network or
%                         archive drive
%     MasconFile ("")     a mascon .nc, exported as SHX_MASCON_FILE and
%                         persisted as preference 'MasconFile'. This
%                         enables the opt-in readMascon check against a
%                         REAL product: the synthetic fixture writes
%                         (lon, lat, time) while the JPL file writes
%                         (time, lat, lon), and a size check passes
%                         either way. Sources: CSR (no login), GSFC, or
%                         JPL PO.DAAC (free Earthdata login)
%     GravisFolder ("")   GravIS aux folder for the validated chains,
%                         exported as SHX_GRAVIS_FOLDER + preference
%                         'GravisFolder' (opt-in real-data tests)
%     DDKFolder ("")      DDK filter binaries folder, exported as
%                         SHX_DDK_FOLDER + preference 'DDKFolder'
%     Quiet (false)       suppress progress output
%     DataFolder ("")  persistent data folder, applied BEFORE any fetcher runs
%     FetchITSG ("none")  "all" additionally downloads every monthly ITSG solution
%     Proxy ("")  per-call proxy URL, e.g. "http://proxy:8080" (empty: MATLAB Web Preferences)
%
%   Output
%     summary struct:
%       root        (1,1) string   toolbox root folder
%       pathAction  (1,1) string   "already-on-path" | "added-temporarily"
%                                  | "added-permanently" | "planned"
%       dataFolder  (1,1) string   download target root
%       plan        (1,:) string   human-readable action list
%       fetched     (1,:) string   files newly downloaded (verified)
%       skipped     (1,:) string   files already present (re-verified)
%       failed      (1,:) string   files that failed (download or parse)
%       testData    (1,:) string   opt-in test-data variables that were
%                                  set, if any (SHX_SERIES_FOLDER,
%                                  SHX_MASCON_FILE). setenv lasts for
%                                  the session only - put the same lines
%                                  in startup.m to keep them
%       ok          (1,1) logical  isempty(failed)
%
%   Examples
%     setup_shAnalysis(Permanent = true, Download = "core")
%     setup_shAnalysis(Download = "filters", DDK = [2 3 5])
%     s = setup_shAnalysis(DryRun = true, Download = "starter")
%     setup_shAnalysis(SeriesFolder = "D:/grace/itsg", ...
%         MasconFile = "D:/grace/GRCTellus_JPL.nc")   % opt-in test data
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-07 (v2.5).
arguments
    opts.Permanent (1,1) logical = false
    opts.Download (1,1) string {mustBeMember(opts.Download, ...
        ["none", "core", "filters", "starter"])} = "none"
    opts.DDK (1,:) double {mustBeInteger, mustBeInRange(opts.DDK, 1, 8)} = 3
    opts.Providers (1,:) string ...
        {mustBeMember(opts.Providers, ["GFZ", "CSR", "JPL"])} = ...
        ["GFZ", "CSR", "JPL"]
    opts.Months (1,:) string = ["2008-04", "2025-12"]
    opts.Nmax (1,1) double {mustBeMember(opts.Nmax, [60, 96])} = 96
    opts.Docs (1,1) logical = false
    opts.DryRun (1,1) logical = false
    opts.DataFolder (1,1) string = ""
    opts.FetchITSG (1,1) string = "none"
    opts.Proxy (1,1) string = ""
    opts.Update (1,1) logical = false
    opts.SeriesFolder (1,1) string = ""
    opts.MasconFile (1,1) string = ""
    opts.GravisFolder (1,1) string = ""
    opts.DDKFolder (1,1) string = ""
    opts.Quiet (1,1) logical = false
end
root = string(fileparts(mfilename('fullpath')));
level = find(opts.Download == ["none", "core", "filters", "starter"]);

% ------------------------------------------------------------- the plan
plan = strings(1, 0);
entries = split(string(path), pathsep);
if ispc
    onPath = any(strcmpi(entries, root));
else
    onPath = any(strcmp(entries, root));
end
if onPath
    plan(end+1) = "path: already on the MATLAB path (" + root + ")";
elseif opts.Permanent
    plan(end+1) = "path: addpath + savepath (permanent): " + root;
else
    plan(end+1) = "path: addpath (this session only): " + root;
end
if level >= 2
    plan(end+1) = "core: shLowLevel.fetchTN -> TN-14 (GSFC) + TN-13 (" + ...
        strjoin(opts.Providers, "/") + "), verified by parse" + ...
        ternary(opts.Update, " (update existing)", "");
end
if level >= 3
    plan(end+1) = "filters: shLowLevel.fetchDDK -> DDK" + ...
        strjoin(string(unique(opts.DDK)), ", DDK") + " (~10 MB each)" + ...
        ternary(opts.Update, " (update existing)", "");
end
if level >= 4
    plan(end+1) = "starter: shLowLevel.fetchITSG -> ITSG n" + opts.Nmax + ...
        " months " + strjoin(opts.Months, ", ") + ...
        ternary(opts.Update, " (update existing)", "");
end
if opts.Docs
    plan(end+1) = "docs: builddocsearchdb(html/) for help-browser search";
end

summary = struct('root', root, 'pathAction', "planned", ...
    'dataFolder', "", 'plan', plan, 'fetched', strings(1, 0), ...
    'skipped', strings(1, 0), 'failed', strings(1, 0), ...
    'testData', strings(1, 0), 'ok', true);

% ---- data folder override BEFORE any fetcher runs (v3.0.0): the
% pre-v3 chicken-and-egg (dataFolder had to be set via a function that
% is only on the path after setup) is gone
if strlength(opts.DataFolder) > 0 && ~opts.DryRun
    shLowLevel.dataFolder(char(opts.DataFolder));
end
if opts.DryRun
    % read-only dataFolder lookup (shLowLevel.dataFolder() would mkdir)
    if ispref('shAnalysis', 'dataFolder')
        summary.dataFolder = string(getpref('shAnalysis', 'dataFolder'));
    else
        summary.dataFolder = string(fullfile(root, 'data'));
    end
    if strlength(opts.SeriesFolder) > 0
        summary.testData(end+1) = "SHX_SERIES_FOLDER = " + opts.SeriesFolder;
    end
    if strlength(opts.MasconFile) > 0
        summary.testData(end+1) = "SHX_MASCON_FILE = " + opts.MasconFile;
    end
    if ~opts.Quiet
        fprintf('setup_shAnalysis DRY RUN - planned actions:\n');
        fprintf('  %s\n', plan);
        if ~isempty(summary.testData)
            fprintf('  would set: %s\n', summary.testData);
        end
    end
    return
end

% ---------------------------------------------------------------- path
if onPath
    summary.pathAction = "already-on-path";
else
    addpath(char(root));
    summary.pathAction = "added-temporarily";
end
if opts.Permanent
    st = savepath;
    if st == 0
        summary.pathAction = "added-permanently";
    else
        warning('shAnalysis:setup:savepathFailed', ...
            ['savepath failed (no write permission on pathdef.m?). ' ...
             'The path is set for this session; to persist it, add\n' ...
             '    addpath(''%s'')\nto your startup.m ' ...
             '(see "doc startup").'], root);
    end
end
if ~opts.Quiet
    fprintf('shAnalysis %s\n', summary.pathAction);
end

% ----------------------------------------------------------- downloads
summary.dataFolder = string(shLowLevel.dataFolder());
if level >= 2
    try
        [~, iTN] = shLowLevel.fetchTN(Providers = opts.Providers, ...
            Update = opts.Update, Proxy = opts.Proxy, Quiet = opts.Quiet);
        summary.fetched = [summary.fetched, iTN.fetched, iTN.updated];
        summary.skipped = [summary.skipped, iTN.skipped];
        summary.failed  = [summary.failed,  iTN.failed];
    catch err
        summary.failed(end+1) = "fetchTN: " + err.message;
    end
end
% ---- optional bulk ITSG download (v3.0.0): FetchITSG = "all" pulls
% every monthly solution (ITSG-Grace2018 + operational, Nmax = 96)
if opts.FetchITSG == "all"
    try
        [~, iI] = shLowLevel.fetchITSG("all", Nmax = 96, Update = opts.Update, ...
            Proxy = opts.Proxy, Quiet = opts.Quiet);
        summary.fetched = [summary.fetched, iI.fetched, iI.updated];
        summary.skipped = [summary.skipped, iI.skipped];
    catch err
        summary.failed(end+1) = "fetchITSG: " + err.message;
    end
elseif opts.FetchITSG ~= "none"
    error('shAnalysis:setup:badFetchITSG', ...
        'FetchITSG must be "none" or "all" (got "%s").', opts.FetchITSG);
end
if level >= 3
    try
        [~, iD] = shLowLevel.fetchDDK(unique(opts.DDK), ...
            Update = opts.Update, Proxy = opts.Proxy, Quiet = opts.Quiet);
        summary.fetched = [summary.fetched, iD.fetched, iD.updated];
        summary.skipped = [summary.skipped, iD.skipped];
    catch err
        summary.failed(end+1) = "fetchDDK: " + err.message;
    end
end
if level >= 4
    try
        [~, iI] = shLowLevel.fetchITSG(opts.Months, Nmax = opts.Nmax, ...
            Update = opts.Update, Proxy = opts.Proxy, Quiet = opts.Quiet);
        summary.fetched = [summary.fetched, string(iI.fetched), ...
            string(iI.updated)];
        summary.skipped = [summary.skipped, string(iI.skipped)];
        if isfield(iI, 'missing') && ~isempty(iI.missing)
            summary.failed = [summary.failed, ...
                "fetchITSG missing: " + strjoin(string(iI.missing), ", ")];
        end
    catch err
        summary.failed(end+1) = "fetchITSG: " + err.message;
    end
end

% ----------------------------------------------- opt-in test data paths
% Two test suites need data far too large to ship as fixtures: the trend
% regression in testScience needs a monthly series, and the mascon
% reader check needs a real .nc. Both read an environment variable and
% filter themselves out when it is empty, so setting them here is the
% one place a user has to think about it.
%
% Data locations: setenv for this session AND setpref('shAnalysis', ...)
% for persistence across sessions - the opt-in tests and chains read
% both (env first). getpref survives restarts; no startup.m needed.
locs = ["SeriesFolder", "SHX_SERIES_FOLDER", "folder"
        "MasconFile",   "SHX_MASCON_FILE",   "file"
        "GravisFolder", "SHX_GRAVIS_FOLDER", "folder"
        "DDKFolder",    "SHX_DDK_FOLDER",    "folder"];
for li = 1:size(locs, 1)
    val = opts.(locs(li, 1));
    if strlength(val) == 0, continue; end
    isOk = (locs(li, 3) == "folder" && isfolder(val)) || ...
           (locs(li, 3) == "file" && isfile(val));
    if ~isOk
        warning("shAnalysis:setup:no" + locs(li, 1), ...
            '%s is not a %s: %s', locs(li, 1), locs(li, 3), val);
    end
    if ~opts.DryRun
        setenv(locs(li, 2), char(val));
        setpref('shAnalysis', char(locs(li, 1)), char(val));
    end
    summary.testData(end+1) = locs(li, 1) + " = " + val + ...
        "   (env " + locs(li, 2) + " + setpref)"; %#ok<AGROW>
end
if ~isempty(summary.testData) && ~opts.Quiet
    fprintf('data locations (session env + persistent preference):\n');
    fprintf('  %s\n', summary.testData);
end

% ---------------------------------------------------------------- docs
if opts.Docs
    try
        builddocsearchdb(char(fullfile(root, 'html')));
        if ~opts.Quiet, fprintf('help-browser search index built.\n'); end
    catch err
        warning('shAnalysis:setup:docsFailed', ...
            'builddocsearchdb failed: %s', err.message);
    end
end

summary.ok = isempty(summary.failed);
if ~opts.Quiet
    fprintf(['setup done: %d fetched, %d already present, %d failed. ' ...
             'Data folder: %s\n'], numel(summary.fetched), ...
        numel(summary.skipped), numel(summary.failed), summary.dataFolder);
    if ~summary.ok
        fprintf('  failures:\n');
        fprintf('    %s\n', summary.failed);
    end
end
end

function s = ternary(tf, a, b)
if tf, s = a; else, s = b; end
end
