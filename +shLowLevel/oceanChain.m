function [out, rep] = oceanChain(folder, opts)
%OCEANCHAIN Global-ocean mass chain: barystatic series from GSM + GravIS.
%
%   [OUT, REP] = shLowLevel.oceanChain(FOLDER, kn = KN, OceanMask = OC)
%   runs the ocean sibling of the validated chains (roadmap item 7):
%   gravisL2B corrections (C20/C30/C21/S21, degree 1, minus NFIL mean,
%   GIA rate ON by default - the ocean-floor GIA correction matters for
%   barystatic trends) -> global EWH synthesis on a GridStep graticule
%   -> area-weighted ocean-mean series [cm] -> trend + annual fit. The
%   residual after the per-pixel trend+seasonal fit is what remains for
%   the noise proxy: this chain SEPARATES the real open-ocean signal
%   (barystatic rise + seasonal circulation) from the residual that
%   oceanRMS should see - the guide-V6 finding ("the trend's open-ocean
%   RMS is mostly real barystatic signal, not error") as an API.
%
%   Declared limitation: the dealiasing product GAD is NOT restored by
%   default (pass GADFolder= to add it where available); AOD1B carries
%   no secular trend by construction, so the TREND is interpretable
%   without GAD while sub-annual ocean variability is incomplete.
%
%   Inputs
%     folder  (1 x 1) string  monthly .gfc folder (GSM-2_*, any centre)
%
%   Options
%     kn (required)      (N x 2 | N x 1) load Love numbers
%     OceanMask (required) (fn handle | nLat x nLon logical) ocean
%                        decision @(lat, lon) -> logical, MUST be false
%                        over land (oceanRMS contract; the chains'
%                        documented 5 Gt/yr latitude-band lesson)
%     gravisFolder ("")  GravIS aux folder; "" = shipped data/gravis
%     GIAFile (default)  "" disables the GIA rate (see gravisL2B)
%     GIAEpoch (2011)    (1 x 1) GIA reference epoch
%     GridStep (1)       (1 x 1) synthesis graticule [deg]
%     SpanEnd ([])       (1 x 1) keep epochs <= SpanEnd [decimal years]
%     GADFolder ("")     folder of GAD-2_*.gfc; when set, epochs with a
%                        begin-matched GAD file get it added back and
%                        REP records the uncovered count
%     Quiet (true)       (1 x 1) suppress progress output
%
%   Outputs
%     out  (1 x 1) struct with fields
%       epochs    (T x 1) decimal years
%       oceanMean (T x 1) area-weighted ocean-mean EWH [cm]
%       trend     (1 x 1) barystatic trend [mm/yr]
%       trendGt   (1 x 1) the same as mass rate [Gt/yr] over the mask
%       trendSigma (1 x 1) OLS 1-sigma of the trend [mm/yr] (no AR
%                 correction - stated, not hidden)
%       ampAnnual (1 x 1) annual amplitude of the ocean mean [mm]
%       sigMon    (1 x 1) open-ocean residual RMS after the per-pixel
%                 trend+seasonal fit [m] (the honest noise proxy)
%       oceanArea (1 x 1) mask area [m^2]
%     rep  (1 x 1) struct  steps, version, dropped/uncovered epochs
%
%   Example
%     kn = readmatrix("loadLoveNumbers_Gegout97.txt", FileType = "text", ...
%         NumHeaderLines = 2);
%     oc = @(la, lo) abs(la) <= 66 & ~yourLandFcn(la, lo);
%     [out, rep] = shLowLevel.oceanChain("E:/series/COSTG", kn = kn, ...
%         OceanMask = oc);
%     plot(out.epochs, out.oceanMean)   % barystatic curve [cm]
%
%   Error identifiers
%     shLowLevel:oceanChain:missingKn  kn not supplied
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.9.0).
arguments
    folder (1,1) string
    opts.kn double = []
    opts.OceanMask = []
    opts.gravisFolder (1,1) string = ""
    opts.GIAFile string {mustBeScalarOrEmpty} = missing
    opts.GIAEpoch (1,1) double = 2011
    opts.GridStep (1,1) double {mustBePositive} = 1
    opts.SpanEnd double {mustBeScalarOrEmpty} = []
    opts.GADFolder (1,1) string = ""
    opts.Quiet (1,1) logical = true
end
if isempty(opts.kn)
    error('shLowLevel:oceanChain:missingKn', ...
        'oceanChain requires load Love numbers: pass kn = ... explicitly.');
end
R = 6378136.3; GM = 3.986004415e14; Re = 6371e3;
% ---- corrected series (shared core)
l2bArgs = {'SpanEnd', opts.SpanEnd, 'GIAEpoch', opts.GIAEpoch};
if ~ismissing(opts.GIAFile), l2bArgs = [l2bArgs, {'GIAFile', opts.GIAFile}]; end
[ts, repL] = shLowLevel.gravisL2B(folder, opts.gravisFolder, l2bArgs{:});
ep = ts.epochs(:); T = numel(ep);
steps = repL.steps;
% ---- optional GAD restore (begin-matched on filenames)
nGad = 0;
if strlength(opts.GADFolder) > 0
    [ts, nGad] = local_addGAD(ts, opts.GADFolder);
    steps(end+1) = sprintf("GAD restored for %d/%d epochs (folder %s)", ...
        nGad, T, opts.GADFolder);
end
% ---- global synthesis + masks
st = opts.GridStep;
lat = (-90+st/2 : st : 90-st/2)'; lon = (st/2 : st : 360-st/2)';
if isa(opts.OceanMask, 'function_handle')
    [LO, LA] = meshgrid(lon, lat);
    mk = logical(opts.OceanMask(LA, LO));
else
    mk = logical(opts.OceanMask);
end
kn = opts.kn; if size(kn, 2) > 1, kn = kn(:, 2); end
nmaxS = size(ts.Cs, 1) - 1;
E = zeros(numel(lat), numel(lon), T); Pl = [];
for k = 1:T
    if isempty(Pl)
        [E(:,:,k), ~, ~, Pl] = shLowLevel.shSynthesis(ts.Cs(:,:,k), ...
            ts.Ss(:,:,k), GM, R, lat, lon, 'quantity','ewh', 'kn',kn, 'nmin',0);
    else
        E(:,:,k) = shLowLevel.shSynthesis(ts.Cs(:,:,k), ts.Ss(:,:,k), ...
            GM, R, lat, lon, 'quantity','ewh', 'kn',kn, 'nmin',0, 'P',Pl);
    end
end
steps(end+1) = sprintf("EWH synthesis %d epochs, %g-deg grid, nmax %d", ...
    T, st, nmaxS);
% ---- ocean-mean series + fits
dphi = deg2rad(st); dlam = deg2rad(st);
[~, LA2] = meshgrid(lon, lat);
Apix = (Re^2) * cosd(LA2) * dphi * dlam;
wA = Apix(mk); oceanArea = sum(wA);
X = reshape(E, [], T);
om = (wA' * X(mk(:), :))' / oceanArea * 100;              % cm
A6 = [ones(T,1), ep-mean(ep), cos(2*pi*ep), sin(2*pi*ep), ...
      cos(4*pi*ep), sin(4*pi*ep)];
x = A6 \ om;
r = om - A6*x;
sig = sqrt(sum(r.^2) / (T - 6));
Sxx = sum((ep-mean(ep)).^2);
trend = x(2) * 10;                                        % cm/yr -> mm/yr
trendSigma = sig / sqrt(Sxx) * 10;
ampAnnual = hypot(x(3), x(4)) * 10;                       % mm
trendGt = (x(2)/100) * oceanArea * 1000 / 1e12;           % m/yr*m^2*rho/1e12
% ---- per-pixel residual RMS (the honest noise proxy after separation)
coef = A6 \ X(mk(:), :)';
residRMS = sqrt(mean((X(mk(:), :)' - A6*coef).^2, 1));
sigMon = sqrt(median(residRMS.^2));                       % robust vs coasts
steps(end+1) = sprintf(...
    "ocean mean: trend %+.2f mm/yr (%+.1f Gt/yr), amp %.1f mm, sigMon %.4f m", ...
    trend, trendGt, ampAnnual, sigMon);
out = struct('epochs', ep, 'oceanMean', om, 'trend', trend, ...
    'trendGt', trendGt, 'trendSigma', trendSigma, 'ampAnnual', ampAnnual, ...
    'sigMon', sigMon, 'oceanArea', oceanArea);
rep = struct('steps', steps(:), 'version', shLowLevel.version(), ...
    'nEpochs', T, 'nGadRestored', nGad);
if ~opts.Quiet, fprintf('%s\n', steps); end
end

function [ts, nGad] = local_addGAD(ts, gadFolder)
% add back GAD-2 monthly means, matched on the begin date in the name
df = dir(fullfile(char(gadFolder), 'GAD-2_*.gfc'));
tok = regexp({df.name}, 'GAD-2_(\d{4})(\d{3})-', 'tokens', 'once');
yd = @(y) 365 + double(mod(y,4)==0 & (mod(y,100)~=0 | mod(y,400)==0));
begG = nan(numel(tok), 1);
for k = 1:numel(tok)
    y = str2double(tok{k}{1}); d = str2double(tok{k}{2});
    begG(k) = y + (d-1)/yd(y);
end
nGad = 0; nmax = size(ts.Cs, 1) - 1;
for k = 1:ts.nEpochs
    [dmin, j] = min(abs(begG - (ts.epochs(k) - 15/365)));  % ~mid - half month
    if isempty(dmin) || dmin > 0.05, continue; end
    g = shCoefficients.read(fullfile(df(j).folder, df(j).name));
    nm = min(nmax, size(g.C, 1) - 1);
    ts.Cs(1:nm+1, 1:nm+1, k) = ts.Cs(1:nm+1, 1:nm+1, k) + g.C(1:nm+1, 1:nm+1);
    ts.Ss(1:nm+1, 1:nm+1, k) = ts.Ss(1:nm+1, 1:nm+1, k) + g.S(1:nm+1, 1:nm+1);
    nGad = nGad + 1;
end
end
