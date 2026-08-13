function [files, info] = fetchITSGBackground(months, opts)
%FETCHITSGBACKGROUND Download ITSG monthly-mean background models.
%
%   FILES = shLowLevel.fetchITSGBackground("2018-06") downloads the
%   monthly means of the background models the ITSG processing reduced
%   from the observations, as .gfc (about 1 MB each): "dealiasing"
%   (atmosphere + ocean - the ITSG counterpart of restoring AOD),
%   "oceanTide", "earthTide", "poleTide" and "oceanPoleTide". Release
%   routing follows fetchITSG: months before 2017-07 from
%   ITSG-Grace2018, from 2018-06 on from ITSG-Grace_operational.
%
%   NOT the AOD1B GAX family: ITSG background models are per-release
%   monthly means of the models ITSG actually used; the AOD1B
%   GAA/GAB/GAC/GAD split (shLowLevel.fetchGAX) does not exist here.
%   For restoring ocean signal under an ITSG series, "dealiasing" is
%   the closest counterpart of GAC (atmosphere + ocean, global) - it
%   is NOT a GAD (ocean-only bottom pressure) substitute.
%
%   Inputs
%     months  (1 x k string | numeric years) "YYYY-MM" strings or
%             year vectors, as in fetchITSG
%
%   Options
%     Products (["dealiasing"]) (1 x k string) product selection,
%                        one subfolder per product. Availability is
%                        release-dependent (verified live 2026-08-12):
%                        BOTH eras carry dealiasing, earthTide,
%                        oceanTide, poleTide, oceanPoleTide; the
%                        GRACE era (ITSG-Grace2018) additionally
%                        splits atmosphere, ocean,
%                        oceanBottomPressure and provides c20,
%                        degree1, glacialIsostaticAdjustment,
%                        hydrology. A product absent in an era 404s
%                        and lands in info.missing - loudly, never
%                        silently
%     Dest ("")          destination folder; "" =
%                        dataFolder/series/itsg/background (v3.16)
%     Release ("")       "" routes by month; a release name pins it
%     Update (false)     re-download files that already exist
%     MaxFiles (Inf)     downloads cap per product
%     BudgetSec (Inf)    wall-clock cap with clean partial stop
%     MaxFailures (5)    consecutive-failure cap, stops loudly
%     PauseSec (1.5)     polite pause between downloads [s]
%     RetryAfterCap (30) single capped retry wait on HTTP 429 [s]
%     Downloader ("websave") "websave" | "httpFetch"
%     BaseURL ("https://ftp.tugraz.at/pub/ITSG/GRACE")  server base
%     Quiet (false)      suppress progress output
%
%   Outputs
%     files (n x 1) string  local .gfc paths (downloaded or present)
%     info  (1 x 1) struct  per-product nListed/nDownloaded/nSkipped/
%           nFailed/nRemaining, missing months, releases used
%
%   Example
%     f = shLowLevel.fetchITSGBackground(["2018-06", "2018-07"]);
%     g = shCoefficients.read(f(1));   % plain gfc, C(n+1, m+1)
%
%   Error identifiers
%     shLowLevel:fetchITSGBackground:badProduct  unknown product name
%     shLowLevel:fetchITSGBackground:noMonths    months resolves empty
%
%   Data courtesy of TU Graz, ftp.tugraz.at/pub/ITSG/GRACE
%   (Mayer-Guerr et al.); cite the ITSG series when publishing.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.12.0).
arguments
    months
    opts.Products (1,:) string = "dealiasing"
    opts.Dest (1,1) string = ""
    opts.Release (1,1) string = ""
    opts.Update (1,1) logical = false
    opts.MaxFiles (1,1) double = Inf
    opts.BudgetSec (1,1) double {mustBePositive} = Inf
    opts.MaxFailures (1,1) double {mustBePositive} = 5
    opts.PauseSec (1,1) double {mustBeNonnegative} = 1.5
    opts.RetryAfterCap (1,1) double {mustBePositive} = 30
    opts.Downloader (1,1) string ...
        {mustBeMember(opts.Downloader, ["websave", "httpFetch"])} = "websave"
    opts.BaseURL (1,1) string = "https://ftp.tugraz.at/pub/ITSG/GRACE"
    opts.Quiet (1,1) logical = false
end
known = ["dealiasing", "oceanTide", "earthTide", "poleTide", ...
    "oceanPoleTide", "atmosphere", "ocean", "oceanBottomPressure", ...
    "c20", "degree1", "glacialIsostaticAdjustment", "hydrology"];
bad = setdiff(opts.Products, known);
if ~isempty(bad)
    error('shLowLevel:fetchITSGBackground:badProduct', ...
        'unknown product(s): %s (use %s).', strjoin(bad, ', '), ...
        strjoin(known, '/'));
end
mm = normalizeMonthList(months, "shLowLevel:fetchITSGBackground");
if isempty(mm)
    error('shLowLevel:fetchITSGBackground:noMonths', ...
        'months resolves to an empty set.');
end
dest = opts.Dest;
if strlength(dest) == 0
    dest = fullfile(shLowLevel.dataFolder(), "series", "itsg", "background");
end
rel = strings(numel(mm), 1);
for k = 1:numel(mm)
    if strlength(opts.Release) > 0
        rel(k) = opts.Release;
    elseif mm(k) < "2017-07"
        rel(k) = "ITSG-Grace2018";
    else
        rel(k) = "ITSG-Grace_operational";
    end
end
files = strings(0, 1);
pinfo = struct('product', {}, 'nListed', {}, 'nDownloaded', {}, ...
    'nSkipped', {}, 'nFailed', {}, 'nRemaining', {});
for p = opts.Products
    urls = strings(numel(mm), 1); names = strings(numel(mm), 1);
    for k = 1:numel(mm)
        names(k) = sprintf("model_%s_%s.gfc", p, mm(k));
        urls(k) = opts.BaseURL + "/" + rel(k) + ...
            "/monthly/monthly_background/" + names(k);
    end
    [fs, st] = fetchFileSet(urls, names, fullfile(char(dest), char(p)), ...
        Update = opts.Update, MaxFiles = opts.MaxFiles, ...
        BudgetSec = opts.BudgetSec, MaxFailures = opts.MaxFailures, ...
        PauseSec = opts.PauseSec, RetryAfterCap = opts.RetryAfterCap, ...
        Downloader = opts.Downloader, Quiet = opts.Quiet);
    files = [files; fs]; %#ok<AGROW>
    pinfo(end+1) = struct('product', p, 'nListed', st.nListed, ...
        'nDownloaded', st.nDownloaded, 'nSkipped', st.nSkipped, ...
        'nFailed', st.nFailed, 'nRemaining', st.nRemaining); %#ok<AGROW>
    if ~opts.Quiet
        fprintf('%s: %d listed, %d downloaded, %d present, %d failed\n', ...
            p, st.nListed, st.nDownloaded, st.nSkipped, st.nFailed);
    end
end
missing = mm(~arrayfun(@(m) any(contains(files, m)), mm));
info = struct('products', pinfo, 'missing', missing, ...
    'releases', unique(rel));
end
