function [out, rep] = twsChain(folder, gravisFolder, opts)
%TWSCHAIN The validated GravIS terrestrial-water-storage basin chain.
%
%   [OUT, REP] = shLowLevel.twsChain(FOLDER, GRAVISFOLDER, kn=KN)
%   reproduces guide chapter V8: gravisL2B corrections WITH the
%   ICE-6G_D GIA rate (the GravIS TWS product carries this correction
%   even though its Technical Note omits it - subtracting it took the
%   trend RMS from 0.225 to 0.032 cm/yr), a DDK3 filter as the declared
%   stand-in for GravIS' VDK5/VDK3 blend, and cos-weighted 1-degree
%   basin means in cm EWH computed by per-basin bounding-box synthesis
%   with Legendre recycling (about 1 ms per epoch and basin). With the
%   tested defaults, the COST-G RL02.1 series and the GravIS river
%   basins, eleven major basins reproduce the GravIS portal series at a
%   median amplitude ratio of 1.001 and 0.032 cm/yr trend RMS over the
%   full 252-month span. Every input is exchangeable.
%
%   Inputs
%     folder        (1,1) string  monthly GSM-2_*.gfc folder
%     gravisFolder  (1,1) string  GravIS aux folder (see gravisL2B);
%                   must also hold BasinFile and, for the default
%                   filter, the DDK binary unless full paths are given
%
%   Options
%     kn         (:,:) double  REQUIRED load Love numbers
%     BasinFile  ("basins_rivbas.json")  GravIS GeoJSON of the target
%                basins ([lon lat] rings, converted internally;
%                https://gravis.gfz.de/basins/rivbas)
%     Basins     ([])  string array of feature names to process; []
%                processes every feature in BasinFile
%     Filter     ("DDK3")  "none" | "gaussN" | "DDKn" (n = 1..8) | a W
%                struct from shLowLevel.readDDK
%     DDKFolder  ("")  folder holding the Wbd_2-120.a_* binaries for
%                the "DDKn" shorthand (default: GRAVISFOLDER, then the
%                toolbox data folder)
%     GIAFile    ("GRAVIS-2B_COSTG_0200_GIA_ICE-6G_D_VM5a_0001.gz")
%                "" disables the GIA step - and reproduces the 7x
%                larger trend RMS documented in V8
%     SpanEnd    (Inf)  the TWS validation uses the full record
%     GridStep   (1)    basin-box synthesis step [deg]
%     GM         (3.986004415e14), R (6378136.3)
%     Quiet      (false)
%
%   Outputs
%     out  (1,1) struct  name (K,1 string), epochs (T,1), series
%          (K x T double, basin means [cm EWH]), trend (K,1 [cm/yr]),
%          amplitude (K,1 [cm], annual), phase (K,1 [rad])
%     rep  (1,1) struct  steps, filter, l2b (gravisL2B report),
%          version, created
%
%   Example
%     kn = readmatrix("loadLoveNumbers_Gegout97.txt", FileType="text", ...
%         NumHeaderLines=2);
%     [out, rep] = shLowLevel.twsChain("E:/series/COSTG", "E:/GravIS", ...
%         kn = kn, Basins = ["Amazonas", "Congo", "Danube"]);
%     plot(out.epochs, out.series); legend(out.name)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.8.8).
arguments
    folder (1,1) string
    gravisFolder (1,1) string
    opts.kn double = []
    opts.BasinFile (1,1) string = "basins_rivbas.json"
    opts.Basins string = strings(0, 1)
    opts.Filter = "DDK3"
    opts.DDKFolder (1,1) string = ""
    opts.GIAFile (1,1) string = "GRAVIS-2B_COSTG_0200_GIA_ICE-6G_D_VM5a_0001.gz"
    opts.SpanEnd (1,1) double = Inf
    opts.GridStep (1,1) double {mustBePositive} = 1
    opts.GM (1,1) double = 3.986004415e14
    opts.R (1,1) double = 6378136.3
    opts.Quiet (1,1) logical = false
end
if isempty(opts.kn)
    error('shLowLevel:twsChain:noKn', ...
        ['Load Love numbers are required (kn=). The chain will not ' ...
         'assume a reference frame for you - see fetchLoveNumbers.']);
end
steps = strings(1, 0);
% ---- corrected series (GIA on by default: guide V8)
[ts, l2b] = shLowLevel.gravisL2B(folder, gravisFolder, ...
    GIAFile = opts.GIAFile, SpanEnd = opts.SpanEnd, Quiet = opts.Quiet);
ep = ts.epochs(:); T = ts.nEpochs; nmax = ts.nmax;
% ---- filter
Cs = ts.Cs; Ss = ts.Ss; filtName = "custom W";
if isstruct(opts.Filter)
    ts = ts.applyDDK(opts.Filter); Cs = ts.Cs; Ss = ts.Ss;
elseif opts.Filter == "none"
    filtName = "none";
elseif startsWith(opts.Filter, "gauss")
    rkm = double(extractAfter(opts.Filter, "gauss"));
    w = shLowLevel.shGaussianWeights(nmax, rkm); w = w(:);
    for k = 1:T
        Cs(:,:,k) = Cs(:,:,k) .* w; Ss(:,:,k) = Ss(:,:,k) .* w;
    end
    filtName = opts.Filter;
elseif startsWith(opts.Filter, "DDK")
    ddkMap = ["1d14p","1d13p","1d12p","5d11p","1d11p","5d10p","1d10p","5d9p"];
    n = double(extractAfter(opts.Filter, "DDK"));
    if ~(n >= 1 && n <= 8)
        error('shLowLevel:twsChain:badFilter', 'DDKn needs n in 1..8, got %s.', opts.Filter);
    end
    fn = "Wbd_2-120.a_" + ddkMap(n) + "_4";
    cand = [fullfile(char(opts.DDKFolder), char(fn)); ...
            fullfile(char(gravisFolder), char(fn)); ...
            fullfile(shLowLevel.dataFolder(), 'ddk', char(fn))];
    fp = "";
    for c = 1:size(cand, 1)
        if isfile(strtrim(cand(c, :))), fp = strtrim(cand(c, :)); break; end
    end
    if strlength(fp) == 0
        error('shLowLevel:twsChain:missingDDK', ...
            ['DDK binary %s not found (DDKFolder, GRAVISFOLDER, data ' ...
             'folder). Get it from github.com/strawpants/GRACE-filter.'], fn);
    end
    W = shLowLevel.readDDK(fp);
    ts = ts.applyDDK(W); Cs = ts.Cs; Ss = ts.Ss;
    filtName = opts.Filter;
else
    error('shLowLevel:twsChain:badFilter', 'Unknown Filter: %s', string(opts.Filter));
end
steps(end+1) = "filter: " + filtName;
% ---- basins
bf = char(opts.BasinFile);
if ~isfile(bf), bf = fullfile(char(gravisFolder), bf); end
if ~isfile(bf)
    error('shLowLevel:twsChain:missingBasins', ...
        ['Basin GeoJSON not found: %s. Fetch it from ' ...
         'https://gravis.gfz.de/basins/rivbas.'], opts.BasinFile);
end
J = jsondecode(fileread(bf));
allN = arrayfun(@(f) string(f.properties.name), J.features);
want = opts.Basins;
if isempty(want), want = allN; end
miss = setdiff(want, allN);
if ~isempty(miss)
    error('shLowLevel:twsChain:unknownBasin', ...
        'Basin(s) not in %s: %s', opts.BasinFile, strjoin(miss, ', '));
end
K = numel(want);
% ---- per-basin box synthesis (Legendre recycled across epochs)
series = zeros(K, T);
for k = 1:K
    g = J.features(find(allN == want(k), 1)).geometry;
    if strcmp(g.type, 'Polygon'), parts = {squeeze(g.coordinates(1,:,:))};
    else, parts = cellfun(@(c) squeeze(c(1,:,:)), g.coordinates, 'uni', 0);
    end
    Pall = vertcat(parts{:});
    la1 = floor(min(Pall(:,2))) - 1; la2 = ceil(max(Pall(:,2))) + 1;
    latB = (la1 + opts.GridStep/2 : opts.GridStep : la2)';
    c0 = atan2d(mean(sind(Pall(:,1))), mean(cosd(Pall(:,1))));
    lw = mod(Pall(:,1) - c0 + 180, 360) - 180;
    lo1 = floor(min(lw)) - 1; lo2 = ceil(max(lw)) + 1;
    off = (lo1 + opts.GridStep/2 : opts.GridStep : lo2)';
    lonB = mod(c0 + off, 360);
    [LOc, LAb] = meshgrid(off, latB);
    mkB = false(size(LAb));
    for q = 1:numel(parts)
        Pq = parts{q};
        pl = mod(Pq(:,1) - c0 + 180, 360) - 180;
        mkB = mkB | inpolygon(LOc, LAb, pl, Pq(:,2));
    end
    if ~any(mkB(:))
        error('shLowLevel:twsChain:emptyBasin', ...
            'Basin %s covers no grid point at GridStep %g.', want(k), opts.GridStep);
    end
    wB = cosd(LAb); P = [];
    for t = 1:T
        if isempty(P)
            [E, ~, ~, P] = shLowLevel.shSynthesis(Cs(:,:,t), Ss(:,:,t), ...
                opts.GM, opts.R, latB, lonB, 'quantity','ewh','kn',opts.kn,'nmin',0);
        else
            E = shLowLevel.shSynthesis(Cs(:,:,t), Ss(:,:,t), opts.GM, ...
                opts.R, latB, lonB, 'quantity','ewh','kn',opts.kn,'nmin',0,'P',P);
        end
        series(k, t) = sum(E(mkB) .* wB(mkB)) / sum(wB(mkB)) * 100;
    end
end
steps(end+1) = sprintf("basin means: %d basins x %d epochs", K, T);
% ---- fits
A = [ones(T,1), ep-mean(ep), cos(2*pi*ep), sin(2*pi*ep), cos(4*pi*ep), sin(4*pi*ep)];
X = (A \ series')';
trend = X(:, 2);
amplitude = hypot(X(:, 3), X(:, 4));
phase = atan2(X(:, 4), X(:, 3));
if ~opts.Quiet
    for k = 1:K
        fprintf('  %-28s trend %+6.2f cm/yr  amp %5.1f cm\n', ...
            want(k), trend(k), amplitude(k));
    end
end
out = struct('name', want(:), 'epochs', ep, 'series', series, ...
    'trend', trend, 'amplitude', amplitude, 'phase', phase);
rep = struct('steps', steps, 'filter', filtName, 'l2b', l2b, ...
    'version', shLowLevel.version(), 'created', string(datetime('now')));
end
