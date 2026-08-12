function audit_degree1
%AUDIT_DEGREE1 A deliberately messier, PHYSICAL synthetic world for
%   estimateDegree1. Shipped claim: "1.5-4% error at 2 mm noise".
%   Mess added over the reference world: real GRACE land structure,
%   coloured noise, omission (truth n=90, obs n=60), flipped coastline
%   pixels, a 2017-2018 gap, unmodelled dynamic-ocean sloshing.
%   Mass-conserving: land load + eustatic passive ocean; degree-1 truth
%   = the n=1 analysis of the total field with the same operator the
%   estimator inverts (convention cancels by construction).
%   Audit script, Claude (Fable 5), 2026-08-12.
rng(42);
GM = 3.986004415e14; R = 6378136.3;
kn = readmatrix('E:\DATAPOOL\GravityField\loveNumbers\loadLoveNumbers_Gegout97.txt', ...
                FileType='text', NumHeaderLines=2);
latE = (-89:2:89)'; lonE = (0:2:358)';
nT = 60; T = 120;
t = 2004 + (0:T-1)'/12;  t(t>=2017 & t<2018.5) = [];  T = numel(t);

% ---- land pattern: the REAL chain trend field, seasonally scaled
L = load('E:\DATAPOOL\GravityField\GravIS\audit_chain_cache_unf.mat', ...
         'trendGrid','lat','lon');
oc = @(la,lo) trueOceanA(la,lo);
[LO, LA] = meshgrid(L.lon, L.lat);
isOceanFine = oc(LA, LO);
landT = L.trendGrid; landT(isOceanFine) = 0;
% seasonal continental hydrology: cm-scale blobs (Amazon, SE Asia, Congo,
% Mississippi, Ob) -> mm-scale geocentre + ~8 mm eustatic cycle
blob = @(la0,lo0,s,a) a*exp(-(((LA-la0)/s).^2 + ((mod(LO-lo0+180,360)-180)/(s/cosd(la0))).^2));
hyd = blob(-5,295,12,0.12) + blob(20,102,10,0.10) + blob(-2,20,10,0.08) ...
    + blob(38,268,9,0.05) + blob(60,72,10,0.06);
hyd(isOceanFine) = 0;
amp = 0.35;
oceanExtra = @(k) 0.005 * (sin(2*pi*t(k)+0.4) * (sind(LA).^2 - 1/3) ...
                         + cos(2*pi*t(k)) * cosd(LA).^2 .* sind(2*LO)) .* isOceanFine;
% (Y20/Y22-artig: orthogonal zu Y10/Y11/S11 -- Fall (a). Der fruehere
%  cos(lat)cos(lon)-Slosh ~ Y11 ist Fall (b): dokumentierter Silent-Bias.)

dphiF = deg2rad(abs(L.lat(2)-L.lat(1))); dlamF = deg2rad(abs(L.lon(2)-L.lon(1)));
AF = cosd(LA) * dphiF * dlamF;
oceA = sum(AF(isOceanFine));
kf90 = shLowLevel.kernelFactors("ewh", 90, GM, R, kn=kn);
n90 = (0:90)';
sigDeg = 2e-3 * (1 + (n90/25).^2) / (1 + (15/25)^2);   % 2 mm EWH at n=15
Cs = zeros(nT+1, nT+1, T); Ss = zeros(nT+1, nT+1, T);
truth.C10 = zeros(T,1); truth.C11 = zeros(T,1); truth.S11 = zeros(T,1);
[LOe, LAe] = meshgrid(lonE, latE);
oceanE_true = oc(LAe, LOe);
OceanModel = zeros(numel(latE), numel(lonE), T);
for k = 1:T
    s = 1 + amp*cos(2*pi*t(k) - 0.5);
    Lk = s*landT + cos(2*pi*t(k) - 0.6) * hyd;
    u = -sum(Lk(:).*AF(:)) / oceA;                     % eustatic closure
    sig = Lk + u*double(isOceanFine) + oceanExtra(k);  % total EWH [m]
    [Ck, Sk] = shLowLevel.shAnalysisGrid(sig, L.lat, L.lon, 90, ...
        quantity="ewh", GM=GM, R=R, kn=kn);
    truth.C10(k) = Ck(2,1); truth.C11(k) = Ck(2,2); truth.S11(k) = Sk(2,2);
    for n = 2:90
        sc = sigDeg(n+1) / (kf90(n+1)*sqrt(2*n+1));
        Ck(n+1,1:n+1) = Ck(n+1,1:n+1) + sc*randn(1,n+1);
        Sk(n+1,2:n+1) = Sk(n+1,2:n+1) + sc*randn(1,n);
    end
    Cs(:,:,k) = Ck(1:nT+1,1:nT+1); Ss(:,:,k) = Sk(1:nT+1,1:nT+1);
    Cs(1:2,1:2,k) = 0; Ss(1:2,1:2,k) = 0;              % GSM: degrees 2+
    OceanModel(:,:,k) = u * double(oceanE_true);       % passive part only
end
ts = shSeries(Cs, Ss=Ss, Epochs=t, GM=GM, R=R);

% ---- ocean mask with flipped coastline pixels
ocean = oceanE_true;
flip = rand(size(ocean)) < 0.02 & ...
       abs(conv2(double(ocean),ones(3)/9,'same') - ocean) > 0.05;
ocean(flip) = ~ocean(flip);
fprintf('flipped coastline pixels: %d (%.1f%% of grid)\n', nnz(flip), 100*nnz(flip)/numel(flip));

[d1, info] = shLowLevel.estimateDegree1(ts, ocean, kn=kn, ...
    OceanModel=OceanModel, LatDeg=latE, LonDeg=lonE);
score("C10", d1.C10, truth.C10, t, R);
score("C11", d1.C11, truth.C11, t, R);
score("S11", d1.S11, truth.S11, t, R);
fprintf('info: oceanFraction=%.3f  hasModel=%d  cond(max)=%.3g  residRMS(max)=%.4g\n', ...
    info.oceanFraction, info.hasModel, max(info.cond), max(info.residualRMS));

fprintf('-- no-OceanModel control (documented bias) --\n');
d1b = shLowLevel.estimateDegree1(ts, ocean, kn=kn, LatDeg=latE, LonDeg=lonE);
score("C10", d1b.C10, truth.C10, t, R);
end

function score(name, est, tru, t, R)
A = [ones(size(t)), t-mean(t), cos(2*pi*t), sin(2*pi*t)];
xe = A\est; xt = A\tru;
ampE = hypot(xe(3),xe(4)); ampT = hypot(xt(3),xt(4));
mm = sqrt(3)*R*1e3;
fprintf(['%s: annual amp %.3f vs %.3f mm (%+.1f%%) | trend %+.4f vs %+.4f ' ...
    'mm/yr | RMS resid %.3f mm\n'], name, ampE*mm, ampT*mm, ...
    100*(ampE-ampT)/ampT, xe(2)*mm, xt(2)*mm, rms(est-tru)*mm);
end

function tf = trueOceanA(la, lo)
lo = mod(lo+180,360)-180;
land = (la>-35 & la<37  & lo>-17 & lo<52)  | (la>5 & la<75 & lo>25 & lo<180) | ...
       (la>15 & la<72 & lo>-168 & lo<-52)  | (la>-56 & la<13 & lo>-82 & lo<-34) | ...
       (la>-45 & la<-10 & lo>112 & lo<154) | (la<-63) | ...
       (la>59 & la<84 & lo>-73 & lo<-12);
tf = ~land;
end
