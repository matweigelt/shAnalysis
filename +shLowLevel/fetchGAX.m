function [files, info] = fetchGAX(dest, opts)
%FETCHGAX Download AOD1B GAX monthly means (GAA/GAB/GAC/GAD) from ICGEM.
%
%   FILES = shLowLevel.fetchGAX(DEST) downloads the monthly-mean
%   dealiasing products needed for full ocean-mass restoration
%   (Chambers & Willis 2010) as ICGEM-converted .gfc files into
%   DEST/<PRODUCT>/. The AOD1B products are identical across processing
%   centres, so the GFZ series pages (which carry all four products
%   from 2002 on) serve every GSM series including COST-G. Listing and
%   download go through the toolbox's Retry-After-aware fetch layer;
%   files already present are skipped unless Update = true.
%
%   Inputs
%     dest    (1 x 1) string  destination folder; product subfolders
%             are created as needed
%
%   Options
%     Products (["GAD", "GAA"])  (1 x k string) any of GAA/GAB/GAC/GAD;
%              GAD restores the model ocean signal, the ocean mean of
%              GAA removes the atmospheric land-ocean mass term
%     Series (["01_GRACE/GFZ/GFZ Release 06", "01_GRACE/GFZ/GFZ Release 06.3 (GFO)"])
%              (1 x s string) ICGEM catalogue paths
%              "group/center/series" fed to shLowLevel.listICGEM
%     Update (false) (1 x 1) re-download files that already exist
%     MaxFiles (Inf) (1 x 1) stop after this many downloads per
%              product - bounded acceptance runs, not a user knob
%     Downloader ("websave") (1 x 1) "websave" (webread network stack,
%              matches listICGEM) or "httpFetch" (matlab.net.http with
%              Retry-After support); the stacks differ in proxy/TLS
%              behaviour, websave is the robust default for ICGEM
%     BudgetSec (Inf) (1 x 1) wall-clock budget; on expiry the fetch
%              stops CLEANLY, keeps what it has, and reports the
%              remainder (info.nRemaining) - never a hung call
%     MaxFailures (5) (1 x 1) consecutive-failure cap per product;
%              on hitting it the product loop stops with a warning
%              instead of failing every remaining file one by one
%     PauseSec (1.5) (1 x 1) polite pause between downloads [s] - the
%              ICGEM server rate-limits bursts (HTTP 429; observed live)
%     RetryAfterCap (30) (1 x 1) on a 429 the download is retried ONCE
%              after this wait [s]; a second 429 counts as a failure
%              (websave does not expose the Retry-After header, so a
%              fixed cap replaces server-driven waits)
%     Quiet (false)  (1 x 1) suppress progress output
%
%   Outputs
%     files (n x 1) string  full paths of all matching local files
%           (downloaded now or already present)
%     info  (1 x 1) struct  nListed, nDownloaded, nSkipped, nFailed,
%           nRemaining (budget cut) per product, series actually used
%
%   Example
%     f = shLowLevel.fetchGAX("E:/DATAPOOL/GravityField/GAX");
%     [out, rep] = shLowLevel.oceanChain(ser, kn = kn, OceanMask = oc, ...
%         GADFolder = "E:/DATAPOOL/GravityField/GAX/GAD", ...
%         GAAFolder = "E:/DATAPOOL/GravityField/GAX/GAA");
%
%   See fetchSINEX and fetchITSGBackground for the ITSG-side products;
%   the ITSG background models are per-release means, NOT the AOD1B
%   GAX split.
%
%   Error identifiers
%     shLowLevel:fetchGAX:badProduct  product not one of GAA/GAB/GAC/GAD
%     shLowLevel:fetchGAX:noFiles     a requested product listed zero
%                                     files on every series page
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.10.0).
arguments
    dest (1,1) string
    opts.Products (1,:) string = ["GAD", "GAA"]
    opts.Series (1,:) string = ["01_GRACE/GFZ/GFZ Release 06", ...
        "01_GRACE/GFZ/GFZ Release 06.3 (GFO)"]
    opts.Update (1,1) logical = false
    opts.MaxFiles (1,1) double = Inf
    opts.Downloader (1,1) string ...
        {mustBeMember(opts.Downloader, ["websave", "httpFetch"])} = "websave"
    opts.BudgetSec (1,1) double {mustBePositive} = Inf
    opts.MaxFailures (1,1) double {mustBePositive} = 5
    opts.PauseSec (1,1) double {mustBeNonnegative} = 1.5
    opts.RetryAfterCap (1,1) double {mustBePositive} = 30
    opts.Quiet (1,1) logical = false
end
bad = setdiff(upper(opts.Products), ["GAA","GAB","GAC","GAD"]);
if ~isempty(bad)
    error('shLowLevel:fetchGAX:badProduct', ...
        'unknown product(s): %s (use GAA/GAB/GAC/GAD).', strjoin(bad, ', '));
end
prods = upper(opts.Products);
files = strings(0, 1);
pinfo = struct('product', {}, 'nListed', {}, 'nDownloaded', {}, ...
    'nSkipped', {}, 'nFailed', {}, 'nRemaining', {});
% one listing per series page, reused for every product
L = cell(1, numel(opts.Series));
for si = 1:numel(opts.Series)
    T = shLowLevel.listICGEM(Type = "temporal", Series = opts.Series(si));
    L{si} = T;
end
for pi = 1:numel(prods)
    p = prods(pi);
    sub = fullfile(char(dest), char(p));
    if ~isfolder(sub), mkdir(sub); end
    Uall = strings(0,1); Nall = strings(0,1);
    for si = 1:numel(opts.Series)
        T = L{si};
        keep = contains(T.url, "GAX_products/" + p + "/");
        Uall = [Uall; T.url(keep)]; Nall = [Nall; T.name(keep)]; %#ok<AGROW>
    end
    nL = numel(Uall);
    [fs, st] = fetchFileSet(Uall, Nall, string(sub), ...
        Update = opts.Update, MaxFiles = opts.MaxFiles, ...
        BudgetSec = opts.BudgetSec, MaxFailures = opts.MaxFailures, ...
        PauseSec = opts.PauseSec, RetryAfterCap = opts.RetryAfterCap, ...
        Downloader = opts.Downloader, Quiet = opts.Quiet);
    files = [files; fs];
    nD = st.nDownloaded; nS = st.nSkipped; nF = st.nFailed; nRem = st.nRemaining;
    if nL == 0
        error('shLowLevel:fetchGAX:noFiles', ...
            '%s: zero files listed on every series page.', p);
    end
    pinfo(end+1) = struct('product', p, 'nListed', nL, 'nDownloaded', nD, ...
        'nSkipped', nS, 'nFailed', nF, 'nRemaining', nRem); %#ok<AGROW>
    if ~opts.Quiet
        fprintf('%s: %d listed, %d downloaded, %d present, %d failed\n', ...
            p, nL, nD, nS, nF);
    end
end
info = struct('products', pinfo, 'series', opts.Series);
end
