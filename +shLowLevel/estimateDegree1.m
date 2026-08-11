function [d1, info] = estimateDegree1(ts, ocean, opts)
%ESTIMATEDEGREE1 Geocentre coefficients from GRACE itself (Swenson 2008).
%
%   [D1, INFO] = shLowLevel.estimateDegree1(TS, OCEAN) estimates the
%   degree-1 coefficients C10, C11, S11 for every epoch of the series TS,
%   using the method of Swenson, Chambers & Wahr (2008).
%
%   GRACE cannot sense degree 1: those coefficients describe the offset
%   between the centre of mass and the centre of figure, the satellites
%   orbit the centre of mass, and the coefficients are set to zero by
%   definition. Every surface-mass estimate needs them, so they have to
%   come from somewhere - usually a published table (see
%   shLowLevel.readTN13). This function derives them from the data
%   instead, which removes the dependency on somebody else having
%   published a degree-1 series for your exact Level-2 product, and
%   guarantees the result is consistent with the solutions you actually
%   used. It is what GravIS and geogravL3 do.
%
%   The argument: a surface mass field is land plus ocean. GRACE observes
%   degrees 2 and up of the total. The degree-1 terms are then whatever,
%   when added, makes the OCEAN part of the field agree with an ocean
%   model - a small least-squares problem per epoch, solved over the
%   ocean domain with area weighting. A constant term is carried
%   alongside the three patterns to absorb any mass imbalance between
%   data and ocean model; without it that offset lands on C10 (a 153%
%   error on the reference problem) because C10 has a large ocean mean.
%
%   Inputs
%     ts         (1,1) shSeries  the solutions, degrees 2+ (a GSM series;
%                any existing degree 1 is IGNORED, not added to)
%     ocean      (nLat x nLon) logical, or a function handle
%                f(latDeg, lonDeg) -> logical. true = ocean. This is
%                user-supplied, like Love numbers: a wrong ocean domain
%                changes the answer silently
%
%   Options
%     kn ([])  load Love numbers, degrees 0..nmax (REQUIRED - the
%             degree-1 value is frame dependent and must not be assumed;
%             see shLowLevel.readLoveNumbers for the frame conversion)
%     OceanModel ([])  (nLat x nLon x T) double  expected ocean mass in
%             EWH [m] per epoch, or (nLat x nLon) for a static field, or
%             [] for none. OMITTING IT BIASES THE RESULT: whatever ocean
%             signal the model would have explained is then attributed
%             to the geocentre instead (validated: up to 177% error on
%             the reference problem)
%     LatDeg (-89:2:89), LonDeg (0:2:358)  the evaluation grid. It only
%             has to resolve the ocean domain and the three smooth
%             degree-1 patterns, so a coarse grid is fine and fast
%     Nmax (NaN)  truncate the observed field before evaluation
%             (NaN: the series' own nmax)
%
%   Outputs
%     d1         (1,1) struct  fields: epoch (T,1 double),
%                C10, C11, S11 (T,1 double) - the same layout
%                shLowLevel.readTN13 returns, so it drops straight into
%                shSeries.addDegree1 / shCoefficients.addDegree1
%     info       (1,1) struct  fields: cond (T,1 double, condition
%                number of the column-equilibrated design matrix per
%                epoch (equilibrated so the number reflects the GEOMETRY
%                of the ocean domain, not the units of the columns) - a rising value
%                warns that the ocean domain is too small or too
%                lopsided BEFORE the noise shows it), oceanFraction
%                (1,1 double, area fraction used), nPixels, hasModel
%                (1,1 logical), residualRMS (T,1 double, of the ocean
%                fit in EWH [m])
%
%   Validated in tools/dev/validate_degree1.py: exact recovery of a known
%   geocentre from a synthetic land/ocean world, 1.5-4% error at 2 mm
%   noise, a demonstrable bias when the ocean model is omitted, and a
%   rising condition number on a degenerate ocean domain.
%
%   Example
%     LN = shLowLevel.readLoveNumbers("loadLoveNumbers.txt");
%     d1 = shLowLevel.estimateDegree1(ts, isOcean, kn = LN.kn);
%     ts = ts.addDegree1(d1);            % same interface as TN-13
%
%   See also shLowLevel.readTN13, shSeries.addDegree1,
%   shLowLevel.oceanRMS.
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.6.0).
arguments
    ts (1,1) shSeries
    ocean
    opts.kn double = []
    opts.OceanModel double = []
    opts.LatDeg (1,:) double = -89:2:89
    opts.LonDeg (1,:) double = 0:2:358
    opts.Nmax (1,1) double = NaN
end
assert(~isempty(opts.kn), 'shLowLevel:estimateDegree1:noLoveNumbers', ...
    ['Love numbers are required: kn = ... . The degree-1 value is ' ...
     'frame dependent (CE/CF/CM) and this toolbox never assumes one. ' ...
     'See shLowLevel.readLoveNumbers for the frame conversion.']);
lat = opts.LatDeg(:).';
lon = opts.LonDeg(:).';
nmax = opts.Nmax;
if ~isfinite(nmax), nmax = ts.nmax; end
[LO, LA] = meshgrid(lon, lat);
if isa(ocean, 'function_handle')
    isOc = logical(ocean(LA, LO));
else
    assert(isequal(size(ocean), size(LA)), ...
        'shLowLevel:estimateDegree1:badMask', ...
        'ocean must be a function handle or a %d x %d mask.', ...
        numel(lat), numel(lon));
    isOc = logical(ocean);
end
assert(any(isOc(:)), 'shLowLevel:estimateDegree1:emptyOcean', ...
    'The ocean mask selects no points.');

T = ts.nEpochs;
mdl = opts.OceanModel;
hasModel = ~isempty(mdl);
if hasModel
    if ismatrix(mdl)
        assert(isequal(size(mdl), size(LA)), ...
            'shLowLevel:estimateDegree1:badModel', ...
            'OceanModel must be nLat x nLon [x T].');
        mdl = repmat(mdl, 1, 1, T);
    else
        assert(isequal(size(mdl), [size(LA), T]), ...
            'shLowLevel:estimateDegree1:badModel', ...
            'OceanModel must be nLat x nLon x %d.', T);
    end
end

% the three degree-1 surface-mass patterns, in EWH metres
B = degree1Basis(lat, lon, opts.kn);
w = cosd(LA);
% A CONSTANT column goes with the three patterns. The ocean residual
% generally contains a degree-0 (total mass) component - the observed
% field is degrees 2+, any mass imbalance between data and ocean model
% shows up as an offset - and without a constant to absorb it, that
% offset is projected onto C10, which has a large ocean-mean value. On
% the reference problem this alone produced a 153% error. The constant
% is a nuisance parameter: estimated, then discarded.
A = [B{1}(isOc), B{2}(isOc), B{3}(isOc), ones(nnz(isOc), 1)] ...
    .* sqrt(w(isOc));
% Column-equilibrate before solving. The constant column has a completely
% different scale from the three basis patterns, and an unscaled cond(A)
% then reports the units rather than the geometry (3.8e7 on a problem
% whose geometry is fine). Scaling is undone on the solution, so the
% answer is unchanged and the diagnostic becomes meaningful.
colScale = vecnorm(A, 2, 1);
colScale(colScale == 0) = 1;
A = A ./ colScale;

C10 = zeros(T, 1); C11 = zeros(T, 1); S11 = zeros(T, 1);
cnd = zeros(T, 1); rms = zeros(T, 1);
for k = 1:T
    g = ts.at(k);
    obs = g.synthesis(lat, lon, quantity = "ewh", kn = opts.kn, ...
        nmin = 2, nmax = nmax);
    if hasModel
        d = mdl(:, :, k) - obs;
    else
        d = -obs;
    end
    b = d(isOc) .* sqrt(w(isOc));
    x = (A \ b) ./ colScale(:);            % undo the equilibration
    C10(k) = x(1); C11(k) = x(2); S11(k) = x(3);   % x(4): mass nuisance
    cnd(k) = cond(A);
    rms(k) = sqrt(mean((A * (x .* colScale(:)) - b).^2));
end
d1 = struct('epoch', ts.epochs(:), 'C10', C10, 'C11', C11, 'S11', S11);
info = struct('cond', cnd, 'oceanFraction', sum(w(isOc)) / sum(w(:)), ...
    'nPixels', nnz(isOc), 'hasModel', hasModel, 'residualRMS', rms);
if ~hasModel
    warning('shLowLevel:estimateDegree1:noOceanModel', ...
        ['No OceanModel given. Whatever real ocean signal a model ' ...
         'would have explained is now attributed to the geocentre ' ...
         'instead, which biases the estimate (up to 177%% on the ' ...
         'reference problem in tools/dev/validate_degree1.py).']);
end
end

% ------------------------------------------------------------- helpers
function B = degree1Basis(lat, lon, kn)
%DEGREE1BASIS The C10, C11, S11 surface-mass patterns in EWH [m].
%   A unit degree-1 Stokes coefficient maps to surface density through
%   the Wahr et al. (1998) kernel at n = 1; the same constants as
%   shLowLevel.kernelFactors, kept explicit here because only n = 1 is
%   needed and the basis must match the synthesis exactly.
R = 6378136.3;
rhoAve = 5517;
rhoWater = 1000;
assert(numel(kn) >= 2, 'shLowLevel:estimateDegree1:shortLoveNumbers', ...
    'kn must cover at least degrees 0 and 1.');
kf = R * rhoAve / (3 * rhoWater) * 3 / (1 + kn(2));
th = deg2rad(90 - lat(:));
lam = deg2rad(lon(:).');
P10 = sqrt(3) * cos(th);
P11 = sqrt(3) * sin(th);
B = {kf * P10 .* ones(size(lam)), ...
     kf * P11 .* cos(lam), ...
     kf * P11 .* sin(lam)};
end
