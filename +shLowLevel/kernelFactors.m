function kernel = kernelFactors(quantity, nmax, GM, R, opts)
%KERNELFACTORS Degree-dependent spectral factors for SH synthesis/analysis.
%
%   KERNEL = shLowLevel.kernelFactors(QUANTITY, NMAX, GM, R, kn ([])=..., ...) returns
%   the (NMAX+1)x1 factor f_n such that a field of QUANTITY is
%       q(lat,lon) = sum_n f_n * sum_m Pbar_nm (C_nm cos + S_nm sin).
%   Shared by shLowLevel.shSynthesis (multiply) and shLowLevel.shAnalysisGrid (divide);
%   the single source of truth for the quantity definitions.
%
%   Quantities (v2.3 additions marked *):
%     'geoid'                 [m]        R
%     'potential'             [m^2/s^2]  GM/R
%     'gravity_anomaly'       [m/s^2]    (GM/R^2)(n-1)
%     'gravity_disturbance'   [m/s^2]    (GM/R^2)(n+1)
%     'gravity_gradient_rr' * [1/s^2]    (GM/R^3)(n+1)(n+2)   (T_rr;
%                             1 Eotvos = 1e-9 1/s^2)
%     'ewh'                   [m]        R rho_ave (5517)/(3 rho_w) (2n+1)/(1+kn)
%     'surface_density'     * [kg/m^2]   R rho_ave/3 (2n+1)/(1+kn)
%                             (= rho_w * ewh kernel; requires kn)
%     'bottom_pressure'     * [Pa]       g0 * surface_density kernel,
%                             g0 = GM/R^2 (ocean-bottom / surface
%                             pressure equivalent; requires kn)
%     'deformation_up'      * [m]        R hn ([])/(1+kn): elastic vertical
%                             load deformation (Wahr et al. 1998; the
%                             GRACE <-> GNSS uplift comparison quantity;
%                             requires kn AND hn). Horizontal components
%                             need the gradient synthesis - see
%                             shLowLevel.shSynthesisDeformation.
%
%   Inputs
%     quantity  char/string  see above
%     nmax      (1,1) double
%     GM, R     (1,1) double [m^3/s^2], [m]
%     opts.kn        (>=nmax+1,1) double  load Love numbers k'_n
%                    (NEVER hardcoded - user-supplied)
%     opts.hn        (>=nmax+1,1) double  load Love numbers h'_n
%                    (deformation_up only; user-supplied)
%     opts.rho_ave   (1,1) double = 5517   [kg/m^3]
%     opts.rho_water (1,1) double = 1000   [kg/m^3]
%     opts.Height    (1,1) double = 0      [m] evaluation height above R:
%                    multiplies the potential-type kernels by the upward-
%                    continuation attenuation (R/r)^p with r = R+Height
%                    (p = n+1 potential/geoid-as-potential-ratio... see
%                    below; p = n+2 anomaly/disturbance; p = n+3 T_rr).
%                    Python-validated against -d/dr and d2/dr2 of the
%                    continued potential (2e-10 / 7e-7 = check's own
%                    truncation). Surface quantities (ewh,
%                    surface_density, bottom_pressure, deformation_up,
%                    geoid) reject Height (shSynthesis:heightInvalid) -
%                    they are surface concepts.
%   Outputs
%     kernel    (nmax+1,1) double  (zeros where the quantity carries no
%               information, e.g. n=1 for gravity_anomaly)
%
%   Claude (Fable 5), 2026-08-07 (v2.3 same day).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    quantity {mustBeTextScalar}
    nmax (1,1) double {mustBeInteger, mustBeNonnegative}
    GM (1,1) double
    R (1,1) double
    opts.kn double = []
    opts.hn double = []
    opts.rho_ave (1,1) double = 5517
    opts.rho_water (1,1) double = 1000
    opts.Height (1,1) double {mustBeNonnegative} = 0
end

% Love-number contract (audit F-13/F-19): a short vector used to die as
% MATLAB:badsubscript, and 1+kn = 0 (the CM-frame k1 = -1 convention,
% present in GROOPS' own ak135 files) produced a silent Inf kernel.
if ~isempty(opts.kn)
    assert(numel(opts.kn) >= nmax + 1, 'shLowLevel:kernelFactors:knTooShort', ...
        'kn has %d entries but degrees 0..%d need %d (kn(1) is degree 0).', ...
        numel(opts.kn), nmax, nmax + 1);
    den = 1 + opts.kn(1:nmax+1);
    kb = find(~isfinite(den) | den == 0, 1);
    assert(isempty(kb), 'shLowLevel:kernelFactors:badLoveNumbers', ...
        ['1+kn = %g at degree %d - the load kernel divides by it. ' ...
         '(A CM-frame set with k1 = -1 cannot be used below degree 2.)'], ...
        den(max(kb,1)), kb - 1);
end
if ~isempty(opts.hn)
    assert(numel(opts.hn) >= nmax + 1, 'shLowLevel:kernelFactors:hnTooShort', ...
        'hn has %d entries but degrees 0..%d need %d.', numel(opts.hn), nmax, nmax + 1);
end

n = (0:nmax)';
needKn = @(q) assert(~isempty(opts.kn), 'shSynthesis:missingLoveNumbers', ...
    ['%s requires degree-dependent load Love numbers k''_n.\n' ...
     'Pass them explicitly via ''kn'', e.g. from a validated ' ...
     'PREM/Wahr(1998) table.'], q);
switch lower(char(quantity))
    case 'geoid'
        kernel = R * ones(nmax+1, 1);
    case 'potential'
        kernel = (GM/R) * ones(nmax+1, 1);
    case 'gravity_anomaly'
        kernel = (GM/R^2) * (n - 1);
    case 'gravity_disturbance'
        kernel = (GM/R^2) * (n + 1);
    case 'gravity_gradient_rr'
        kernel = (GM/R^3) * (n + 1) .* (n + 2);
    case 'ewh'
        needKn('EWH');
        kernel = (R * opts.rho_ave) / (3 * opts.rho_water) ...
            * (2*n + 1) ./ (1 + opts.kn(1:nmax+1));
    case 'surface_density'
        needKn('surface_density');
        kernel = (R * opts.rho_ave) / 3 ...
            * (2*n + 1) ./ (1 + opts.kn(1:nmax+1));
    case 'bottom_pressure'
        needKn('bottom_pressure');
        kernel = (GM/R^2) * (R * opts.rho_ave) / 3 ...
            * (2*n + 1) ./ (1 + opts.kn(1:nmax+1));
    case 'deformation_up'
        needKn('deformation_up');
        assert(~isempty(opts.hn), 'shSynthesis:missingLoveNumbers', ...
            ['deformation_up requires load Love numbers h''_n in ' ...
             'addition to k''_n. Pass them explicitly via ''hn''.']);
        kernel = R * opts.hn(1:nmax+1) ./ (1 + opts.kn(1:nmax+1));
    otherwise
        error('shSynthesis:badQuantity', 'Unknown quantity: %s', quantity);
end

if opts.Height > 0
    q = lower(char(quantity));
    r = R + opts.Height;
    switch q
        case 'potential'
            kernel = kernel .* (R/r).^(n + 1);
        case {'gravity_anomaly', 'gravity_disturbance'}
            kernel = kernel .* (R/r).^(n + 2);
        case 'gravity_gradient_rr'
            kernel = kernel .* (R/r).^(n + 3);
        otherwise
            error('shSynthesis:heightInvalid', ...
                ['Height applies to potential-type quantities only ' ...
                 '(potential, gravity_anomaly, gravity_disturbance, ' ...
                 'gravity_gradient_rr); ''%s'' is a surface concept.'], q);
    end
end
end
