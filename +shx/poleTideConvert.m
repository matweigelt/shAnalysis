function varargout = poleTideConvert(varargin)
%POLETIDECONVERT Convert C21/S21 between mean-pole (pole tide) conventions.
%
%   G2 = shx.poleTideConvert(G, From="IERS2010", To="IERS2018") adjusts
%   the C21/S21 coefficients of an shCoefficients object (epoch taken
%   from G) from one mean-pole convention to another. Processing centers
%   remove the (solid + ocean) pole tide using a MEAN POLE MODEL, and
%   the IERS 2010 conventional model (cubic until 2010.0, linear after)
%   differs from the 2018 secular-pole update (linear, adopted for
%   GRACE/GRACE-FO RL06) by up to ~1e-10 in C21/S21 with a strong TREND
%   component - silently mixing conventions between solutions biases
%   C21/S21 trend studies at the mm-EWH level. This function makes the
%   conversion explicit.
%
%   [C2, S2] = shx.poleTideConvert(C, S, EPOCH, ...) is the raw form on
%   coefficient triangles (only the (3,2) = C21/S21 entries change).
%
%   Mean-pole models [mas], dt = t - 2000:
%     IERS2010  t<2010: xm = 55.974+1.8243 dt+0.18413 dt^2+0.007024 dt^3
%                       ym = 346.346+1.7896 dt-0.10729 dt^2-0.000908 dt^3
%               t>=2010: xm = 23.513+7.6141 dt ; ym = 358.891-0.6287 dt
%     IERS2018  xm = 55.0+1.677 dt ; ym = 320.5+3.460 dt (secular pole)
%   Adjustment (published solution A -> convention B), with
%   dm1 = xmB - xmA, dm2 = ymA - ymB [arcsec]:
%     solid (IERS):  dC21 = -1.333e-9 (dm1 + 0.0115 dm2)
%                    dS21 = -1.333e-9 (dm2 - 0.0115 dm1)
%     ocean (Desai): dC21 = -2.1778e-10 (dm1 - 0.01724 dm2)
%                    dS21 = -2.1778e-10 (dm2 + 0.03365 dm1)
%   Python-validated: A->A = 0 and A->B = -(B->A) exactly; the dS21
%   trend matches the dominant-term prediction to < 1%; magnitudes
%   1e-11..1e-10 over 2005-2020.
%
%   Options
%     From (1,1) string = "IERS2010"   convention of the input solution
%     To (1,1) string = "IERS2018"     target convention
%     Mode (1,1) string = "both"       "solid" | "ocean" | "both"
%     SolidCoef (1,1) double = -1.333e-9      overridable constants
%     OceanCoef (1,1) double = -2.1778e-10
%
%   Outputs
%     G2         (1,1) shCoefficients  converted object (object input;
%                                      history appended)
%     C2, S2     (nmax+1 x nmax+1) double  converted triangles (raw input)
%
%   Example
%     g18 = shx.poleTideConvert(g, From = "IERS2010", To = "IERS2018");
%     dC21 = g18.C(3, 2) - g.C(3, 2)     % ~1e-10 level, trend-like in t
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-10 (v2.7.0).
if nargin >= 1 && isa(varargin{1}, 'shCoefficients')
    g = varargin{1};
    if ~isfinite(g.epoch)
        error('shx:poleTideConvert:noEpoch', ...
            'Object has no finite epoch; use the (C, S, EPOCH) form.');
    end
    [C2, S2] = doConvert(g.C, g.S, g.epoch, varargin{2:end});
    % mutate through the class API (properties are SetAccess = private);
    % setCoefficient also appends the history entry
    varargout{1} = g.setCoefficient(2, 1, C2(3, 2), S2(3, 2));
else
    [varargout{1:2}] = doConvert(varargin{:});
end
end

function [C2, S2] = doConvert(C, S, epoch, opts)
arguments
    C double
    S double
    epoch (1,1) double {mustBeFinite}
    opts.From (1,1) string {mustBeMember(opts.From, ...
        ["IERS2010", "IERS2018"])} = "IERS2010"
    opts.To (1,1) string {mustBeMember(opts.To, ...
        ["IERS2010", "IERS2018"])} = "IERS2018"
    opts.Mode (1,1) string {mustBeMember(opts.Mode, ...
        ["solid", "ocean", "both"])} = "both"
    opts.SolidCoef (1,1) double = -1.333e-9
    opts.OceanCoef (1,1) double = -2.1778e-10
end
if size(C, 1) < 3 || ~isequal(size(C), size(S))
    error('shx:poleTideConvert:badSize', ...
        'C and S must be equal-size triangles with nmax >= 2.');
end
[xmA, ymA] = meanPole(opts.From, epoch);
[xmB, ymB] = meanPole(opts.To, epoch);
dm1 = xmB - xmA;                          % m1A - m1B  [arcsec]
dm2 = ymA - ymB;                          % m2A - m2B  [arcsec]
dC = 0; dS = 0;
if opts.Mode == "solid" || opts.Mode == "both"
    dC = dC + opts.SolidCoef * (dm1 + 0.0115 * dm2);
    dS = dS + opts.SolidCoef * (dm2 - 0.0115 * dm1);
end
if opts.Mode == "ocean" || opts.Mode == "both"
    dC = dC + opts.OceanCoef * (dm1 - 0.01724 * dm2);
    dS = dS + opts.OceanCoef * (dm2 + 0.03365 * dm1);
end
C2 = C; S2 = S;
C2(3, 2) = C(3, 2) + dC;
S2(3, 2) = S(3, 2) + dS;
end

function [xm, ym] = meanPole(model, t)
% mean/secular pole [arcsec]
dt = t - 2000.0;
switch model
    case "IERS2010"
        if t < 2010.0
            xm = 55.974 + 1.8243*dt + 0.18413*dt^2 + 0.007024*dt^3;
            ym = 346.346 + 1.7896*dt - 0.10729*dt^2 - 0.000908*dt^3;
        else
            xm = 23.513 + 7.6141*dt;
            ym = 358.891 - 0.6287*dt;
        end
    case "IERS2018"
        xm = 55.0 + 1.677*dt;
        ym = 320.5 + 3.460*dt;
end
xm = xm * 1e-3; ym = ym * 1e-3;           % mas -> arcsec
end
