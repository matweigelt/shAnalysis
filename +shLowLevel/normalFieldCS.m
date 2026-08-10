function [Cn, info] = normalFieldCS(nmax, opts)
%NORMALFIELDCS Even zonal harmonics of an ellipsoidal normal field.
%
%   [CN, INFO] = shLowLevel.normalFieldCS(NMAX, System="WGS84") returns the
%   fully normalized even zonal coefficients Cbar_{2n,0} of the normal
%   gravity field, COMPUTED from the defining constants (GM (NaN), a (NaN), f (NaN),
%   omega (NaN)) via the closed-form theory (Heiskanen & Moritz 1967,
%   sect. 2-9) - no coefficient tables anywhere:
%
%     e^2 = f(2-f),  e' = e/sqrt(1-e^2),  m = w^2 a^2 b / GM
%     q0  = 1/2 [ (1 + 3/e'^2) atan(e') - 3/e' ]
%     J2  = (e^2/3) (1 - (2/15) m e'/q0)
%     J2n = (-1)^{n+1} (3 e^{2n})/((2n+1)(2n+3)) (1 - n + 5n J2/e^2)
%     Cbar_{2n,0} = -J2n / sqrt(4n+1)
%
%   Python-validated against the published NIMA TR8350.2 WGS84 values:
%   J2 = 1.082629821313e-3 and Cbar20/40/60/80 to all published digits.
%
%   Inputs
%     nmax    (1,1) double   maximum degree of the output
%   Options
%     System ("WGS84")  "WGS84" | "GRS80": presets for the DEFINING
%                       constants below (citable definitions, not data)
%     GM, a, f, omega   override any defining constant. GRS80 defines
%                       J2 rather than f; the preset carries the
%                       standard derived 1/f = 298.257222101 (Moritz
%                       2000) so one code path serves both.
%   Outputs
%     Cn    (nmax+1,1) double  Cbar_{n,0} (zero at odd/omitted degrees;
%           Cn(1) = 1 represents the central GM/r term - subtract it
%           ONLY together with a matching degree-0 of the model, see
%           subtractNormalField)
%     info  struct: GM, a, f, omega, J2, e2, m, q0, system
%
%   The normal field is w.r.t. ITS OWN (GM, a): rescale before mixing
%   with model coefficients (shLowLevel.rescaleGMR does; g.subtractNormalField
%   handles the whole chain). Permanent tide: the ellipsoid is a
%   geometric convention (tide-free by construction); the C20 of your
%   model carries its own tide system (ICGEM header 'tide_system') -
%   difference ~4.2e-9 in Cbar20 between zero-tide and tide-free. Not
%   converted silently; check your model's convention.
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    nmax (1,1) double {mustBeInteger, mustBeNonnegative}
    opts.System (1,1) string ...
        {mustBeMember(opts.System, ["WGS84","GRS80"])} = "WGS84"
    opts.GM (1,1) double = NaN
    opts.a (1,1) double = NaN
    opts.f (1,1) double = NaN
    opts.omega (1,1) double = NaN
end
switch opts.System
    case "WGS84"   % NIMA TR8350.2 defining constants
        d = struct('GM', 3.986004418e14, 'a', 6378137.0, ...
            'f', 1/298.257223563, 'omega', 7.292115e-5);
    case "GRS80"   % Moritz (2000); f is the standard derived constant
        d = struct('GM', 3.986005e14, 'a', 6378137.0, ...
            'f', 1/298.257222101, 'omega', 7.292115e-5);
end
for fn = ["GM", "a", "f", "omega"]
    if ~isnan(opts.(fn)), d.(fn) = opts.(fn); end
end
b  = d.a * (1 - d.f);
e2 = d.f * (2 - d.f);
e  = sqrt(e2);
ep = e / sqrt(1 - e2);
m  = d.omega^2 * d.a^2 * b / d.GM;
q0 = 0.5 * ((1 + 3/ep^2) * atan(ep) - 3/ep);
J2 = (e2/3) * (1 - (2/15) * (m * ep / q0));
Cn = zeros(nmax + 1, 1);
Cn(1) = 1;                                  % central term GM/r
for k = 1:floor(nmax/2)
    J2n = (-1)^(k+1) * (3 * e^(2*k)) / ((2*k+1)*(2*k+3)) ...
        * (1 - k + 5*k*J2/e2);
    Cn(2*k + 1) = -J2n / sqrt(4*k + 1);
end
info = struct('GM', d.GM, 'a', d.a, 'f', d.f, 'omega', d.omega, ...
    'J2', J2, 'e2', e2, 'm', m, 'q0', q0, 'system', opts.System);
end
