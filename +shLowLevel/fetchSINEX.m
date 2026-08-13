function [files, info] = fetchSINEX(months, opts)
%FETCHSINEX Download ITSG monthly normal-equation SINEX from TU Graz.
%
%   FILES = shLowLevel.fetchSINEX("2018-06") downloads the monthly
%   normal-equation/solution SINEX of the ITSG series - the only
%   public per-month SINEX source (COST-G combines on the
%   normal-equation level internally but does NOT distribute per-month
%   SINEX or covariances). Release routing follows fetchITSG: months
%   before 2017-07 come from ITSG-Grace2018, months from 2018-06 on
%   from ITSG-Grace_operational. Feed the files to
%   shLowLevel.readSINEX (verified against exactly this product).
%
%   SIZE WARNING, read before fetching: ONE monthly n96 SINEX is about
%   460 MB gzipped (verified live 2026-08-12); the full 257-month
%   series is on the order of 120 GB. MONTHS is therefore a required
%   argument - there is deliberately no "all" convenience here, unlike
%   fetchITSG where a month is 1 MB.
%
%   Inputs
%     months  (1 x k string | numeric years) "YYYY-MM" strings or
%             year vectors, as in fetchITSG
%
%   Options
%     Dest ("")          destination folder; "" = dataFolder/itsg_sinex
%     Nmax (96)          (1 x 1) 96 or 120 - the two server variants
%     Release ("")       "" routes by month (fetchITSG rule);
%                        "ITSG-Grace2018" or "ITSG-Grace_operational"
%                        pin one release
%     Update (false)     re-download files that already exist
%     MaxFiles (Inf)     downloads cap - bounded acceptance runs
%     BudgetSec (Inf)    wall-clock cap with clean partial stop
%     MaxFailures (5)    consecutive-failure cap, stops loudly
%     PauseSec (1.5)     polite pause between downloads [s]
%     RetryAfterCap (30) single capped retry wait on HTTP 429 [s]
%     Downloader ("websave") "websave" | "httpFetch" (Retry-After-aware,
%                        with its own websave transport fallback)
%     BaseURL ("https://ftp.tugraz.at/pub/ITSG/GRACE")  server base
%     Quiet (false)      suppress progress output
%
%   Outputs
%     files (n x 1) string  local .snx.gz paths (downloaded or present)
%     info  (1 x 1) struct  nListed, nDownloaded, nSkipped, nFailed,
%           nRemaining, missing (requested months absent on the
%           server), releases used
%
%   Example
%     f = shLowLevel.fetchSINEX(["2018-06", "2018-07"], Nmax = 96);
%     snx = shLowLevel.readSINEX(f(1), Only = "estimate");
%
%   Error identifiers
%     shLowLevel:fetchSINEX:badNmax    Nmax not 96 or 120
%     shLowLevel:fetchSINEX:noMonths   months resolves to empty
%
%   Data courtesy of TU Graz, ftp.tugraz.at/pub/ITSG/GRACE
%   (Mayer-Guerr et al.); cite the ITSG series when publishing.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.12.0).
arguments
    months
    opts.Dest (1,1) string = ""
    opts.Nmax (1,1) double = 96
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
if ~any(opts.Nmax == [96, 120])
    error('shLowLevel:fetchSINEX:badNmax', ...
        'Nmax must be 96 or 120 (the two server variants), got %g.', ...
        opts.Nmax);
end
mm = normalizeMonthList(months, "shLowLevel:fetchSINEX");
if isempty(mm)
    error('shLowLevel:fetchSINEX:noMonths', ...
        'months resolves to an empty set - one 460 MB file per month is fetched deliberately, never "all".');
end
dest = opts.Dest;
if strlength(dest) == 0
    dest = fullfile(shLowLevel.dataFolder(), "itsg_sinex");
end
% release routing (the fetchITSG rule)
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
urls = strings(numel(mm), 1); names = strings(numel(mm), 1);
for k = 1:numel(mm)
    stem = rel(k);
    if stem == "ITSG-Grace2018", nameStem = "ITSG-Grace2018";
    else, nameStem = "ITSG-Grace_operational"; end
    names(k) = sprintf("%s_n%d_%s.snx.gz", nameStem, opts.Nmax, mm(k));
    urls(k) = opts.BaseURL + "/" + stem + "/monthly/normals_SINEX/" + ...
        sprintf("monthly_n%d/", opts.Nmax) + names(k);
end
[files, st] = fetchFileSet(urls, names, dest, ...
    Update = opts.Update, MaxFiles = opts.MaxFiles, ...
    BudgetSec = opts.BudgetSec, MaxFailures = opts.MaxFailures, ...
    PauseSec = opts.PauseSec, RetryAfterCap = opts.RetryAfterCap, ...
    Downloader = opts.Downloader, Quiet = opts.Quiet);
% months without a local file after the run are missing on the server
% (mission gap 2017-07..2018-05, dropouts) or were cut by the caps
missing = mm(~arrayfun(@(m) any(contains(files, m)), mm));
info = struct('nListed', st.nListed, 'nDownloaded', st.nDownloaded, ...
    'nSkipped', st.nSkipped, 'nFailed', st.nFailed, ...
    'nRemaining', st.nRemaining, 'missing', missing, ...
    'releases', unique(rel));
if ~opts.Quiet
    fprintf('SINEX: %d listed, %d downloaded, %d present, %d failed, %d missing\n', ...
        st.nListed, st.nDownloaded, st.nSkipped, st.nFailed, numel(missing));
end
end
