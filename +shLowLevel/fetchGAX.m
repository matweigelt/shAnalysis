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
%     Quiet (false)  (1 x 1) suppress progress output
%
%   Outputs
%     files (n x 1) string  full paths of all matching local files
%           (downloaded now or already present)
%     info  (1 x 1) struct  nListed, nDownloaded, nSkipped, nFailed
%           per product (struct arrays), series actually used
%
%   Example
%     f = shLowLevel.fetchGAX("E:/DATAPOOL/GravityField/GAX");
%     [out, rep] = shLowLevel.oceanChain(ser, kn = kn, OceanMask = oc, ...
%         GADFolder = "E:/DATAPOOL/GravityField/GAX/GAD", ...
%         GAAFolder = "E:/DATAPOOL/GravityField/GAX/GAA");
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
    'nSkipped', {}, 'nFailed', {});
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
    nL = 0; nD = 0; nS = 0; nF = 0;
    for si = 1:numel(opts.Series)
        T = L{si};
        keep = contains(T.url, "GAX_products/" + p + "/");
        U = T.url(keep); N = T.name(keep);
        nL = nL + numel(U);
        for k = 1:numel(U)
            out = fullfile(sub, char(N(k)));
            if isfile(out) && ~opts.Update
                nS = nS + 1; files(end+1, 1) = string(out); %#ok<AGROW>
                continue
            end
            if nD >= opts.MaxFiles, break, end
            try
                shLowLevel.httpFetch(U(k), out);
                nD = nD + 1; files(end+1, 1) = string(out); %#ok<AGROW>
                if ~opts.Quiet && mod(nD, 25) == 0
                    fprintf('  %s: %d files\n', p, nD);
                end
            catch ME
                nF = nF + 1;
                if ~opts.Quiet
                    fprintf('  [fail] %s: %s\n', N(k), ME.identifier);
                end
            end
        end
    end
    if nL == 0
        error('shLowLevel:fetchGAX:noFiles', ...
            '%s: zero files listed on every series page.', p);
    end
    pinfo(end+1) = struct('product', p, 'nListed', nL, 'nDownloaded', nD, ...
        'nSkipped', nS, 'nFailed', nF); %#ok<AGROW>
    if ~opts.Quiet
        fprintf('%s: %d listed, %d downloaded, %d present, %d failed\n', ...
            p, nL, nD, nS, nF);
    end
end
info = struct('products', pinfo, 'series', opts.Series);
end
