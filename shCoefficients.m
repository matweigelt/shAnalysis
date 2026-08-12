classdef shCoefficients
%SHCOEFFICIENTS Spherical harmonic (Stokes) coefficient set - single epoch.
%
%   The single point of access for one gravity field: a GRACE/GRACE-FO
%   monthly solution (GSM), a dealiasing product (GAA/GAB/GAC/GAD), a
%   static model, or any derived field (difference, filtered, climatology
%   component). Immutable value class: every operation returns a NEW
%   object and appends to .history (full provenance).
%
%   Construction
%     g  = shCoefficients.read('GSM-2_2024032-2024060_....gfc');  % gfc/gfct/.gz
%     g  = shCoefficients(C, S, GM=3.986e14, R=6.378e6, Epoch=2024.12);
%     g  = shCoefficients.fromVec(x, idx, like);
%
%   Typical chain
%     d  = gB - gA;                              % monthly change
%     d2 = d.destripe(minOrder=6).gaussian(300); % fluent filtering
%     [grid, lat, lon] = d2.synthesis(-90:90, 0:359, quantity="ewh", kn=kn);
%     spec = d.degreeRMS;  n0 = d.crossover;  d.spectrum;  d.triangle;
%
%   GAX handling: productType is parsed from the L2 filename (GSM/GAA/
%   GAB/GAC/GAD); g + gad restores the background product, with epoch
%   matching enforced.
%
%   Properties (read-only)
%     C, S           (nmax+1)x(nmax+1) double  fully normalized coefficients
%     sigmaC, sigmaS same layout or []          formal errors
%     GM             (1,1) double  [m^3/s^2]
%     R              (1,1) double  [m]
%     nmax           (1,1) double  (dependent, from size(C))
%     epoch          (1,1) double  decimal years, NaN if unknown/mixed
%     productType    (1,1) string  "GSM"|"GAA".."GAD"|"static"|"difference"|...
%     tideSystem     (1,1) string  from the file header, "unknown" otherwise
%     name           (1,1) string
%     header         struct        raw header key:value pairs
%     variableTerms  struct array  gfct time-variable terms (see evalAt)
%     history        string column provenance log
%
%   See also shSeries, shClimatology.
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

properties (SetAccess = private)
    C double
    S double
    sigmaC double = []
    sigmaS double = []
    GM (1,1) double = 3.986004415e14
    R (1,1) double = 6378136.3
    epoch (1,1) double = NaN
    productType (1,1) string = "unknown"
    tideSystem (1,1) string = "unknown"
    name (1,1) string = ""
    header struct = struct()
    variableTerms = []
    history string = string.empty(0,1)
end

properties (Dependent)
    nmax
end

methods
    function obj = set.history(obj, h)
        % normalize to a string column: scalar-seeded (end+1) growth
        % otherwise produces rows and breaks vertcat chains (v2.2 fix)
        h = string(h);                       % chars -> scalar string first
        obj.history = h(:);
    end

    function obj = shCoefficients(C, S, opts)
        %SHCOEFFICIENTS Construct from coefficient matrices.
        %
        %   Inputs
        %     S  (nmax+1 x nmax+1) double  sine coefficients, same layout (S(:,1) = 0)
        %     C  (nmax+1 x nmax+1) double  cosine coefficients, C(n+1, m+1) layout
        %   Outputs
        %     obj        (1 x 1) shCoefficients   immutable coefficient set with GM (3.986004415e14), R (6378136.3), epoch, sigmas, history
        %
        %   Options
        %     SigmaC ([])  formal 1-sigma of the cosine coefficients, same
        %         layout as C ([]: none attached)
        %     SigmaS ([])  formal 1-sigma of the sine coefficients
        %     Epoch (NaN)  decimal year of the solution; required by every
        %         epoch-matched operation (applyTN14, addDegree1, restore)
        %     ProductType ("unknown")  provider product code, e.g. "GSM",
        %         "GAD", "GAA" - arithmetic uses it to catch mixing levels
        %     TideSystem ("unknown")  "tide_free" | "zero_tide" |
        %         "mean_tide" | "unknown"; carried through and written back
        %         out, never converted silently
        %     Name ("")  model label used in displays and plot titles
        %     Header (struct())  parsed gfc header fields, kept verbatim so
        %         writeGFC can round-trip them
        %     VariableTerms ([])  ICGEM 2.0 gfct terms (trend/periodic per
        %         coefficient) evaluated by evalAt
        %     History (string.empty(0,1))  initial processing history; every
        %         operation appends one line
        %
        %   Outputs
        %     obj  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            C double = 0
            S double = 0
            opts.SigmaC double = []
            opts.SigmaS double = []
            opts.GM (1,1) double {mustBePositive} = 3.986004415e14
            opts.R (1,1) double {mustBePositive} = 6378136.3
            opts.Epoch (1,1) double = NaN
            opts.ProductType (1,1) string = "unknown"
            opts.TideSystem (1,1) string = "unknown"
            opts.Name (1,1) string = ""
            opts.Header struct = struct()
            opts.VariableTerms = []
            opts.History string = string.empty(0,1)
        end
        if ~isequal(size(C), size(S)) || size(C,1) ~= size(C,2)
            error('shCoefficients:badInput', ...
                'C and S must be square matrices of equal size.');
        end
        if ~isempty(opts.SigmaC) && ~isequal(size(opts.SigmaC), size(C))
            error('shCoefficients:badInput', 'sigmaC must match size(C).');
        end
        if ~isempty(opts.SigmaS) && ~isequal(size(opts.SigmaS), size(S))
            error('shCoefficients:badInput', 'sigmaS must match size(S).');
        end
        obj.C = C; obj.S = S;
        obj.sigmaC = opts.SigmaC; obj.sigmaS = opts.SigmaS;
        obj.GM = opts.GM; obj.R = opts.R;
        obj.epoch = opts.Epoch;
        obj.productType = opts.ProductType;
        obj.tideSystem = opts.TideSystem;
        obj.name = opts.Name;
        obj.header = opts.Header;
        obj.variableTerms = opts.VariableTerms;
        obj.history = opts.History;
        if isempty(obj.history)
            obj.history = sprintf("created (nmax=%d)", obj.nmax);
        end
    end

    function n = get.nmax(obj)
        n = size(obj.C, 1) - 1;
    end

    % ---------------------------------------------------------- arithmetic
    function out = plus(a, b)
        %PLUS Coefficient addition (e.g. GSM + GAD background restore).
        %
        %   Inputs
        %     b  (1,1) shCoefficients  right operand (same nmax and GM/R)
        %     a  (1,1) shCoefficients  left operand
        %   Outputs
        %     out        (1 x 1) shCoefficients   sum; sigmas RSS; epochs must match
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        [a, b] = shCoefficients.checkPair(a, b, 'plus');
        shCoefficients.checkEpochs(a, b, 0.02, 'plus');
        out = a.combine(b, +1, ...
            shCoefficients.combinedType(a.productType, b.productType, '+'));
    end

    function out = minus(a, b)
        %MINUS Coefficient difference (e.g. monthly change gB - gA).
        %
        %   Inputs
        %     b  (1,1) shCoefficients  right operand (same nmax and GM/R)
        %     a  (1,1) shCoefficients  left operand
        %   Outputs
        %     out        (1 x 1) shCoefficients   difference; sigmas RSS
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        [a, b] = shCoefficients.checkPair(a, b, 'minus');
        out = a.combine(b, -1, "difference");
        out.epoch = NaN;
        out.history(end) = sprintf("difference: (%s|%.4f) - (%s|%.4f)", ...
            a.name, a.epoch, b.name, b.epoch);
    end

    function out = uminus(a)
        %UMINUS Unary minus: negate all coefficients.
        %
        %   Inputs
        %     a  (1,1) shCoefficients  left operand
        %   Outputs
        %     out        (1 x 1) shCoefficients   negated coefficients (sigmas kept)
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        out = a; out.C = -a.C; out.S = -a.S;
        out.history(end+1) = "negated";
    end

    function out = mtimes(a, b)
        %MTIMES Scalar scaling.
        %
        %   Inputs
        %     b  (1,1) double or shCoefficients  scalar factor (either side) - fields cannot be multiplied
        %   Outputs
        %     out        (1 x 1) shCoefficients   scalar-scaled; sigmas scale by |a|
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        if isa(a, 'shCoefficients') && isnumeric(b) && isscalar(b)
            out = a; f = b;
        elseif isa(b, 'shCoefficients') && isnumeric(a) && isscalar(a)
            out = b; f = a;
        else
            error('shCoefficients:badInput', ...
                'mtimes supports scalar * shCoefficients only.');
        end
        out.C = out.C * f; out.S = out.S * f;
        if ~isempty(out.sigmaC), out.sigmaC = out.sigmaC * abs(f); end
        if ~isempty(out.sigmaS), out.sigmaS = out.sigmaS * abs(f); end
        out.history(end+1) = sprintf("scaled by %.6g", f);
    end

    function [rep, h] = compare(obj, other, varargin)
        %COMPARE Full comparison against another solution (v2.6.0).
        %   REP = g.compare(G2, ...) delegates to shLowLevel.compareSolutions
        %   with OBJ as the first solution; all options pass through
        %   (Plot=true for the 4-panel figure).
        %
        %   Inputs
        %     other  (1,1) shCoefficients  solution to compare against OBJ
        %   Outputs
        %     rep        (1 x 1) struct   shLowLevel.compareSolutions report
        %     h          (1 x 1) graphics handle   figure (Plot = true)
        %
        %   Outputs
        %     rep  (1,1) struct  shLowLevel.compareSolutions report
        %     h    (1,1) graphics handle  figure (Plot = true only)
        [rep, h] = shLowLevel.compareSolutions(obj, other, varargin{:});
    end

    function out = times(a, b)
        %TIMES Elementwise scaling .* (scalar or triangle matrix).
        %   W .* G scales every coefficient by the matching entry of the
        %   (nmax+1 x nmax+1) matrix W (per-coefficient weighting or
        %   tapering); a scalar factor delegates to mtimes. Sigmas scale
        %   by |factor| elementwise (v2.5.1).
        %
        %   Inputs
        %     b  (1,1) shCoefficients or double array  elementwise factor
        %   Outputs
        %     out        (1 x 1) shCoefficients   elementwise-scaled coefficients
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        if isa(a, 'shCoefficients') && isnumeric(b)
            out = a; f = b;
        elseif isa(b, 'shCoefficients') && isnumeric(a)
            out = b; f = a;
        else
            error('shCoefficients:badInput', ...
                'times supports numeric .* shCoefficients only.');
        end
        if isscalar(f)
            out = mtimes(out, f);
            return
        end
        if ~isequal(size(f), size(out.C))
            error('shCoefficients:badInput', ...
                'times factor must be scalar or %d x %d.', ...
                size(out.C, 1), size(out.C, 2));
        end
        out.C = f .* out.C; out.S = f .* out.S;
        if ~isempty(out.sigmaC), out.sigmaC = abs(f) .* out.sigmaC; end
        if ~isempty(out.sigmaS), out.sigmaS = abs(f) .* out.sigmaS; end
        out.history(end+1) = "elementwise scaled";
    end

    % ------------------------------------------------------------ editing
    function out = truncate(obj, nmaxNew)
        %TRUNCATE Reduce the maximum degree.
        %
        %   Inputs
        %     nmaxNew  (1,1) double  new maximum degree (must not exceed the current nmax)
        %   Outputs
        %     out        (1 x 1) shCoefficients   truncated to the new nmax
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            nmaxNew (1,1) double {mustBeInteger, mustBeNonnegative}
        end
        if nmaxNew > obj.nmax
            error('shCoefficients:badTruncation', ...
                'Requested nmax %d exceeds available %d.', nmaxNew, obj.nmax);
        end
        out = obj;
        k = nmaxNew + 1;
        out.C = obj.C(1:k, 1:k); out.S = obj.S(1:k, 1:k);
        if ~isempty(obj.sigmaC), out.sigmaC = obj.sigmaC(1:k, 1:k); end
        if ~isempty(obj.sigmaS), out.sigmaS = obj.sigmaS(1:k, 1:k); end
        out.history(end+1) = sprintf("truncated to nmax=%d", nmaxNew);
    end

    function out = setCoefficient(obj, n, m, Cval, Sval, opts)
        %SETCOEFFICIENT Replace a single (n,m) coefficient pair.
        %
        %   Inputs
        %     Sval  (1,1) double  new S(n+1, m+1) value (NaN keeps the current one)
        %     Cval  (1,1) double  new C(n+1, m+1) value
        %   Outputs
        %     out        (1 x 1) shCoefficients   copy with C/S(n+1, m+1) replaced
        %
        %   Options
        %     SigmaC (NaN)  formal 1-sigma for the coefficient being set
        %         (NaN leaves the existing sigma untouched)
        %     SigmaS (NaN)  same for the sine coefficient
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            n (1,1) double {mustBeInteger, mustBeNonnegative}
            m (1,1) double {mustBeInteger, mustBeNonnegative}
            Cval (1,1) double
            Sval (1,1) double = NaN
            opts.SigmaC (1,1) double = NaN
            opts.SigmaS (1,1) double = NaN
        end
        if n > obj.nmax || m > n
            error('shCoefficients:badIndex', ...
                'Invalid (n,m) = (%d,%d) for nmax=%d.', n, m, obj.nmax);
        end
        out = obj;
        out.C(n+1, m+1) = Cval;
        if ~isnan(Sval), out.S(n+1, m+1) = Sval; end
        if ~isnan(opts.SigmaC) || ~isnan(opts.SigmaS)
            % writing either sigma initializes BOTH stacks full-size:
            % sigmas are paired throughout the toolbox and a partially
            % sigma'd object would break shSeries stacking (v2.5.1)
            if isempty(out.sigmaC), out.sigmaC = nan(obj.nmax + 1); end
            if isempty(out.sigmaS), out.sigmaS = nan(obj.nmax + 1); end
        end
        if ~isnan(opts.SigmaC), out.sigmaC(n+1, m+1) = opts.SigmaC; end
        if ~isnan(opts.SigmaS), out.sigmaS(n+1, m+1) = opts.SigmaS; end
        out.history(end+1) = sprintf("set C/S(%d,%d)", n, m);
    end

    function out = applyTN14(obj, tn, opts)
        %APPLYTN14 Replace C20 (and C30 where provided) from a TN-14 table.
        %   OUT = obj.applyTN14(TN) with TN from shLowLevel.readTN14 (or a
        %   filename). Requires a finite obj.epoch; the nearest TN epoch
        %   within opts.Tolerance (default 0.05 yr) is used.
        %   opts.ReplaceC30: "auto" (default: replace when the table has a
        %   non-NaN C30 for that month), "never", "always".
        %   Outputs
        %     out        (1 x 1) shCoefficients   C20 (and C30 when available) replaced by the SLR values, sigmas updated
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            tn
            opts.Tolerance (1,1) double = 0.05
            opts.ReplaceC30 (1,1) string ...
                {mustBeMember(opts.ReplaceC30, ["auto","never","always"])} = "auto"
        end
        if ischar(tn) || isstring(tn), tn = shLowLevel.readTN14(tn); end
        if ~isfinite(obj.epoch)
            error('shCoefficients:noEpoch', ...
                'applyTN14 requires a finite epoch (set via read/constructor).');
        end
        [dmin, k] = min(abs(tn.epoch - obj.epoch));
        if dmin > opts.Tolerance
            error('shCoefficients:epochNotInTable', ...
                'No TN-14 entry within %.3f yr of epoch %.4f.', ...
                opts.Tolerance, obj.epoch);
        end
        out = obj.setCoefficient(2, 0, tn.C20(k), NaN, SigmaC = tn.sigmaC20(k));
        did = "C20";
        wantC30 = (opts.ReplaceC30 == "always") || ...
                  (opts.ReplaceC30 == "auto" && ~isnan(tn.C30(k)));
        if wantC30 && obj.nmax >= 3
            if isnan(tn.C30(k))
                error('shCoefficients:noC30', ...
                    'ReplaceC30="always" but TN-14 has no C30 for this month.');
            end
            out = out.setCoefficient(3, 0, tn.C30(k), NaN, SigmaC = tn.sigmaC30(k));
            did = did + "+C30";
        end
        out.history(end) = sprintf("TN-14 %s replaced (table epoch %.4f)", ...
            did, tn.epoch(k));
    end

    function out = addDegree1(obj, tn, opts)
        %ADDDEGREE1 Insert TN-13 geocenter (degree-1) coefficients.
        %   OUT = obj.addDegree1(TN) with TN from shLowLevel.readTN13 (or a
        %   filename) sets C10, C11, S11 (with their sigmas) from the
        %   TN-13 record nearest to obj.epoch (within opts.Tolerance (0.05),
        %   default 0.05 yr). Without degree 1, EWH/mass grids are
        %   systematically biased - this completes the standard
        %   GSM + TN-14 + TN-13 processing chain.
        %
        %   Inputs   tn   struct from shLowLevel.readTN13, or filename
        %            Tolerance (1,1) double = 0.05 [yr]
        %   Outputs  out  shCoefficients with degree 1 set
        %   Outputs
        %     out        (1 x 1) shCoefficients   degree-1 row set from the TN-13 record nearest to obj.epoch
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            tn
            opts.Tolerance (1,1) double = 0.05
        end
        if ischar(tn) || isstring(tn), tn = shLowLevel.readTN13(tn); end
        if ~isfinite(obj.epoch)
            error('shCoefficients:noEpoch', ...
                'addDegree1 requires a finite epoch (set via read/constructor).');
        end
        if obj.nmax < 1
            error('shCoefficients:badIndex', ...
                'Model has nmax < 1; cannot hold degree-1 coefficients.');
        end
        [dmin, k] = min(abs(tn.epoch - obj.epoch));
        if dmin > opts.Tolerance
            error('shCoefficients:epochNotInTable', ...
                'No TN-13 entry within %.3f yr of epoch %.4f.', ...
                opts.Tolerance, obj.epoch);
        end
        out = obj.setCoefficient(1, 0, tn.C10(k), NaN, SigmaC = tn.sigC10(k));
        out = out.setCoefficient(1, 1, tn.C11(k), tn.S11(k), ...
            SigmaC = tn.sigC11(k), SigmaS = tn.sigS11(k));
        out.history(end) = sprintf( ...
            "TN-13 degree-1 inserted (table epoch %.4f)", tn.epoch(k));
    end

    % ---------------------------------------------------------- filtering
    function out = destripe(obj, opts)
        %DESTRIPE Swenson & Wahr (2006) / P3M6 decorrelation.
        %   OUT = obj.destripe(minOrder=6, polyOrder=3, windowLength=[])
        %   Formal errors are passed through UNCHANGED (the destriping
        %   operator is data-dependent; rigorous error propagation is not
        %   defined for it - use the tvANS filter on an shSeries for a
        %   filter with proper covariance semantics).
        %   Outputs
        %     out        (1 x 1) shCoefficients   Swenson-Wahr decorrelated
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            opts.minOrder (1,1) double = 6
            opts.polyOrder (1,1) double = 3
            opts.windowLength double = []
        end
        out = obj;
        [out.C, out.S] = shLowLevel.shDestripe(obj.C, obj.S, ...
            'minOrder', opts.minOrder, 'polyOrder', opts.polyOrder, ...
            'windowLength', opts.windowLength);
        if isempty(opts.windowLength)
            out.history(end+1) = sprintf("destriped (minOrder=%d, polyOrder=%d)", ...
                opts.minOrder, opts.polyOrder);
        else
            out.history(end+1) = sprintf(...
                "destriped (minOrder=%d, polyOrder=%d, window=%d)", ...
                opts.minOrder, opts.polyOrder, opts.windowLength);
        end
    end

    function out = gaussian(obj, radiusKm)
        %GAUSSIAN Jekeli (1981) isotropic smoothing; sigmas scaled by Wn.
        %
        %   Inputs
        %     radiusKm  (1,1) double  Gaussian filter half-response radius [km]
        %   Outputs
        %     out        (1 x 1) shCoefficients   Gaussian-smoothed; sigmas scaled by Wn
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            radiusKm (1,1) double {mustBePositive}
        end
        out = obj;
        [out.C, out.S, Wn] = shLowLevel.shGaussianFilter(obj.C, obj.S, radiusKm);
        if ~isempty(obj.sigmaC)
            out.sigmaC = obj.sigmaC .* Wn(:);
        end
        if ~isempty(obj.sigmaS)
            out.sigmaS = obj.sigmaS .* Wn(:);
        end
        out.history(end+1) = sprintf("Gaussian smoothed (%g km)", radiusKm);
    end

    % ----------------------------------------------------------- spectral
    function spec = degreeRMS(obj, opts)
        %DEGREERMS Degree variance/RMS/amplitude spectrum (+ error spectrum).
        %   Outputs
        %     spec       struct: n, amp/rms/var, err   degree spectrum incl. formal-error curve
        %
        %   Options
        %     n0 (0)  first degree included in the spectrum
        arguments
            obj
            opts.n0 (1,1) double = 0
        end
        args = {'R', obj.R, 'n0', opts.n0};
        if ~isempty(obj.sigmaC) && ~isempty(obj.sigmaS)
            args = [args, {'sigmaC', obj.sigmaC, 'sigmaS', obj.sigmaS}];
        end
        spec = shLowLevel.shDegreeRMS(obj.C, obj.S, args{:});
    end

    function [n0, nInterp] = crossover(obj, opts)
        %CROSSOVER Signal-vs-formal-error crossover degree.
        %   Outputs
        %     n0 (2)       (1,1) double  first degree with err >= signal (NaN if none)
        %     nInterp  (1,1) double  interpolated crossing degree
        %
        %   Options
        %     n0 (2)  first degree searched for the signal/error crossover
        arguments
            obj
            opts.n0 (1,1) double = 2
        end
        if isempty(obj.sigmaC) || isempty(obj.sigmaS)
            error('shCoefficients:noSigmas', ...
                'crossover requires formal errors (sigmaC/sigmaS).');
        end
        [n0, nInterp] = shLowLevel.shSpectralCrossover(obj.degreeRMS(n0 = opts.n0));
    end

    function h = spectrum(obj, varargin)
        %SPECTRUM Plot the degree-amplitude spectrum. h = obj.spectrum(...)
        %   Outputs
        %     h          (1 x 1) graphics handle   spectrum plot
        %
        %   Outputs
        %     spec  (1,1) struct  shLowLevel.shDegreeRMS output (degree, degRMS,
        %           degAmplitude, cum*, err* when sigmas are present)
        h = shLowLevel.plotSHSpectrum(obj.degreeRMS, varargin{:});
    end

    function h = triangle(obj, varargin)
        %TRIANGLE Plot the degree/order coefficient triangle. h = obj.triangle(...)
        %   Outputs
        %     h          (1 x 1) graphics handle   coefficient triangle plot
        %
        %   Outputs
        %     h  (1,1) graphics handle  axes of the coefficient triangle plot
        h = shLowLevel.plotSHCoeffTriangle(obj.C, obj.S, varargin{:});
    end

    % ------------------------------------------------------------ spatial
    function [grid, lat, lon] = synthesis(obj, latVec, lonVec, opts)
        %SYNTHESIS Spatial synthesis (geoid, gravity, potential, EWH).
        %   [GRID, LAT, LON] = obj.synthesis(LATVEC, LONVEC, quantity ("geoid")=...,
        %   kn ([])=..., nmin (0)=..., nmax (NaN)=..., UseCache (true)=true).
        %   The Legendre functions are served from a verified process-wide
        %   cache (shLowLevel.legendreCached), so monthly time series on a fixed
        %   grid pay the recursion only once - no manual 'P' passing.
        %   EWH requires explicitly supplied load Love numbers kn (never
        %   hardcoded).
        %   Outputs
        %     grid       (nlat x nlon) double   synthesized field
        %     lat        (1 x nlat) double   geocentric latitudes
        %     lon        (1 x nlon) double   longitudes
        %
        %   Options
        %     hn ([])  vertical-deformation Love numbers, degrees 0..nmax (user-supplied)
        %     Height (0)  evaluation height above the reference sphere [m];
        %         applies the upward continuation (R/r)^n
        %     rho_ave (5517)  mean Earth density [kg/m^3] in the EWH kernel
        %     rho_water (1000)  water density [kg/m^3] in the EWH kernel
        %     Method ("auto")  "auto" | "direct" | "fft": FFT along
        %         longitude is much faster on uniform full-circle longitude
        %         vectors; "auto" uses it whenever it is applicable
        %     LatType ("geocentric")  "geocentric" | "geodetic": how the
        %         latitude inputs are interpreted; geodetic values are
        %         converted with Flattening before evaluation
        %     Flattening (1/298.257223563)  flattening of the reference
        %         ellipsoid for the geodetic/geocentric conversion (WGS84;
        %         overridable, never silently assumed)
        %     MaxMemGB (4)  memory budget [GB]; above it the synthesis
        %         streams in latitude bands instead of allocating the full
        %         design matrix (this is what makes nmax 2190 tractable)
        arguments
            obj
            latVec (1,:) double
            lonVec (1,:) double
            opts.quantity (1,1) string = "geoid"
            opts.kn double = []
            opts.hn double = []
            opts.Height (1,1) double = 0
            opts.nmin (1,1) double = 0
            opts.nmax (1,1) double = NaN
            opts.rho_ave (1,1) double = 5517
            opts.rho_water (1,1) double = 1000
            opts.UseCache (1,1) logical = true
            opts.Method (1,1) string ...
                {mustBeMember(opts.Method, ["auto","direct","fft"])} = "auto"
            opts.LatType (1,1) string ...
                {mustBeMember(opts.LatType, ["geocentric","geodetic"])} = "geocentric"
            opts.Flattening (1,1) double = 1/298.257223563
            opts.MaxMemGB (1,1) double {mustBePositive} = 4
        end
        if opts.LatType == "geodetic"
            % map grids are usually geodetic; SH math needs geocentric
            latVec = shLowLevel.geodetic2geocentric(latVec, ...
                Flattening = opts.Flattening);
        end
        nmaxEff = opts.nmax;
        if isnan(nmaxEff), nmaxEff = obj.nmax; end
        args = {'quantity', char(opts.quantity), 'nmin', opts.nmin, ...
            'nmax', nmaxEff, 'rho_ave', opts.rho_ave, ...
            'rho_water', opts.rho_water, 'method', char(opts.Method), ...
            'MaxMemGB', opts.MaxMemGB, 'Height', opts.Height};
        if ~isempty(opts.kn), args = [args, {'kn', opts.kn}]; end
        if ~isempty(opts.hn), args = [args, {'hn', opts.hn}]; end
        % cache only when the full Legendre stack fits the budget;
        % otherwise let shLowLevel.shSynthesis stream latitude bands (v2.2)
        fitsMem = (nmaxEff+1)^2 * numel(latVec) * 8 <= opts.MaxMemGB * 2^30;
        if opts.UseCache && fitsMem
            P = shLowLevel.legendreCached(nmaxEff, deg2rad(latVec(:)'));
            args = [args, {'P', P}];
        end
        [grid, lat, lon] = shLowLevel.shSynthesis(obj.C, obj.S, obj.GM, obj.R, ...
            latVec, lonVec, args{:});
    end

    function out = toReference(obj, opts)
        %TOREFERENCE Convert to another (GM (NaN), R (NaN)) reference (exact).
        %   OUT = g.toReference(GM=3.986004418e14, R=6378137.0)
        %   re-expresses the coefficients via shLowLevel.rescaleGMR; the
        %   physical field is invariant (Python-validated). Sigmas
        %   rescale identically; GM/R properties are updated.
        %   Outputs
        %     out        (1 x 1) shCoefficients   rescaled to the given GM/R
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            opts.GM (1,1) double = NaN
            opts.R (1,1) double = NaN
        end
        GM2 = opts.GM; if isnan(GM2), GM2 = obj.GM; end
        R2 = opts.R;   if isnan(R2),  R2 = obj.R;  end
        out = obj;
        [out.C, out.S, sg] = shLowLevel.rescaleGMR(obj.C, obj.S, ...
            obj.GM, obj.R, GM2, R2, obj.sigmaC, obj.sigmaS);
        out.sigmaC = sg.C; out.sigmaS = sg.S;
        out.GM = GM2; out.R = R2;
        out.history(end+1) = sprintf( ...
            "rescaled to GM=%.10e, R=%.4f", GM2, R2);
    end

    function out = subtractNormalField(obj, opts)
        %SUBTRACTNORMALFIELD Remove the ellipsoidal normal field.
        %   OUT = g.subtractNormalField()            % WGS84
        %   OUT = g.subtractNormalField(System ("WGS84")="GRS80")
        %   Computes the even zonals from the DEFINING constants
        %   (shLowLevel.normalFieldCS, Heiskanen-Moritz closed form; validated
        %   against NIMA TR8350.2 to all published digits), rescales
        %   them from the ellipsoid's own (GM (NaN), a (NaN)) to the object's
        %   (GM, R) - the WGS84 GM and a DIFFER from the ICGEM values,
        %   skipping this step costs millimetres - and subtracts them
        %   including the degree-0 GM ratio term. The result is the
        %   disturbing field T: 'geoid' synthesis then yields proper
        %   +-100 m undulations (Bruns).
        %
        %   Permanent tide: the ellipsoid is tide-free by construction;
        %   your model's C20 carries its own convention (ICGEM header
        %   'tide_system', ~4.2e-9 zero-tide vs tide-free). Not
        %   converted silently - noted in history only.
        %   Outputs
        %     out        (1 x 1) shCoefficients   disturbing field (even zonals of the normal field removed)
        %
        %   Options
        %     f (NaN)  flattening of the normal field (NaN: the System
        %         default; give GM, a, f and omega together to define your own)
        %     omega (NaN)  angular velocity [rad/s] of the normal field
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            opts.System (1,1) string ...
                {mustBeMember(opts.System, ["WGS84","GRS80"])} = "WGS84"
            opts.GM (1,1) double = NaN
            opts.a (1,1) double = NaN
            opts.f (1,1) double = NaN
            opts.omega (1,1) double = NaN
        end
        [CnEll, infoN] = shLowLevel.normalFieldCS(obj.nmax, ...
            System = opts.System, GM = opts.GM, a = opts.a, ...
            f = opts.f, omega = opts.omega);
        % normal field w.r.t. (GM_ell, a_ell) -> object's (GM, R)
        Cell = zeros(obj.nmax + 1); Cell(:, 1) = CnEll;
        Cn = shLowLevel.rescaleGMR(Cell, zeros(obj.nmax + 1), ...
            infoN.GM, infoN.a, obj.GM, obj.R);
        out = obj;
        out.C(:, 1) = obj.C(:, 1) - Cn(:, 1);
        out.history(end+1) = sprintf( ...
            ['normal field subtracted (%s from defining constants, ' ...
             'rescaled GM %.6e->%.6e, a %.1f->%.1f; permanent-tide ' ...
             'convention of the model unchanged)'], opts.System, ...
            infoN.GM, obj.GM, infoN.a, obj.R);
    end

    function h = map(obj, latDeg, lonDeg, opts)
        %MAP Synthesize and plot in one call (shLowLevel.plotSHMap).
        %   H = g.map(-89:89, 0:359, quantity ("geoid")="ewh", kn ([])=kn) - all
        %   synthesis options plus the plotSHMap display options.
        %
        %   Inputs
        %     lonDeg  (1,:) double  longitudes [deg]
        %     latDeg  (1,:) double  latitudes [deg, geocentric]
        %   Outputs
        %     h          (1 x 1) graphics handle   synthesized map plot
        %
        %   Options
        %     hn ([])  vertical-deformation Love numbers, degrees 0..nmax (user-supplied)
        %     Height (0)  evaluation height above the reference sphere [m]
        %     nmin (0)  lowest degree included in the map
        %     Projection ("plate")  "plate" (plate carree) | "hammer"
        %     Coast (true)  draw coastlines
        %     CLim ([])  color limits [lo hi]; [] uses a robust symmetric
        %         scale from the data percentiles
        %     Units ("")  colorbar label
        %     Title ("")  axes title
        %     ax ([])  target axes handle; [] creates a new figure
        arguments
            obj
            latDeg (1,:) double = -89:2:89
            lonDeg (1,:) double = 0:2:358
            opts.quantity (1,1) string = "geoid"
            opts.kn double = []
            opts.hn double = []
            opts.Height (1,1) double = 0
            opts.nmin (1,1) double = 0
            opts.Projection (1,1) string = "plate"
            opts.Coast (1,1) logical = true
            opts.CLim double = []
            opts.Units (1,1) string = ""
            opts.Title (1,1) string = ""
            opts.ax = []
        end
        grid = obj.synthesis(latDeg, lonDeg, quantity = opts.quantity, ...
            kn = opts.kn, hn = opts.hn, Height = opts.Height, ...
            nmin = opts.nmin);
        ttl = opts.Title;
        if strlength(ttl) == 0
            ttl = sprintf("%s (%s)", obj.name, opts.quantity);
        end
        h = shLowLevel.plotSHMap(grid, latDeg, lonDeg, ...
            Projection = opts.Projection, Coast = opts.Coast, ...
            CLim = opts.CLim, Units = opts.Units, Title = ttl, ...
            ax = opts.ax);
    end

    function out = fan(obj, rDegKm, rOrdKm)
        %FAN Han fan filter (degree x order Gaussian) - shLowLevel.shFanFilter.
        %
        %   Inputs
        %     rOrdKm  (1,1) double  fan filter order-direction radius [km]
        %     rDegKm  (1,1) double  fan filter degree-direction radius [km]
        %   Outputs
        %     out        (1 x 1) shCoefficients   fan-filtered
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            rDegKm (1,1) double {mustBePositive}
            rOrdKm (1,1) double {mustBePositive}
        end
        out = obj;
        [out.C, out.S] = shLowLevel.shFanFilter(obj.C, obj.S, rDegKm, ...
            rOrdKm, R = obj.R);
        out.history(end+1) = sprintf( ...
            "fan filtered (deg %g km, ord %g km)", rDegKm, rOrdKm);
    end

    function [up, north, east] = deformation(obj, latDeg, lonDeg, opts)
        %DEFORMATION Elastic load deformation up/north/east [m].
        %   [UP, NORTH, EAST] = g.deformation(LAT, LON, kn=, hn=, ln=,
        %       Mode ("grid")="grid"|"points", nmin (1)=1, LatType ("geocentric")="geocentric")
        %   from the object's (residual!) coefficients - remove a mean
        %   field first. Love numbers are always user-supplied; see
        %   shLowLevel.shSynthesisDeformation for formulas and validation.
        %
        %   Inputs
        %     lonDeg  (1,:) double  longitudes [deg]
        %     latDeg  (1,:) double  latitudes [deg, geocentric]
        %   Outputs
        %     up         (nlat x nlon | npts) double   vertical deformation [m]
        %     north      same size   north component [m]
        %     east       same size   east component [m]
        %
        %   Options
        %     Flattening (1/298.257223563)  ellipsoid flattening for the
        %         geodetic/geocentric conversion of the station latitudes
        arguments
            obj
            latDeg (1,:) double
            lonDeg (1,:) double
            opts.kn double
            opts.hn double
            opts.ln double
            opts.nmin (1,1) double = 1
            opts.Mode (1,1) string ...
                {mustBeMember(opts.Mode, ["grid","points"])} = "grid"
            opts.LatType (1,1) string ...
                {mustBeMember(opts.LatType, ["geocentric","geodetic"])} = "geocentric"
            opts.Flattening (1,1) double = 1/298.257223563
        end
        if opts.LatType == "geodetic"
            latDeg = shLowLevel.geodetic2geocentric(latDeg, ...
                Flattening = opts.Flattening);
        end
        [up, north, east] = shLowLevel.shSynthesisDeformation(obj.C, obj.S, ...
            obj.R, latDeg, lonDeg, kn = opts.kn, hn = opts.hn, ...
            ln = opts.ln, nmin = opts.nmin, Mode = opts.Mode);
    end

    function out = applyDDK(obj, W)
        %APPLYDDK Apply a DDK (order-block-diagonal) decorrelation filter.
        %   OUT = obj.applyDDK(W) with W from shLowLevel.readDDK (struct, .mat,
        %   or documented ASCII exchange format). Degrees not covered by
        %   the filter pass through unchanged; sigmas are NOT propagated
        %   through the filter (set to NaN - use shLowLevel.mcPropagate for
        %   filtered uncertainties).
        %   Outputs
        %     out        (1 x 1) shCoefficients   DDK-decorrelated
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            W
        end
        W = shLowLevel.readDDK(W);
        out = obj;
        [out.C, out.S] = shLowLevel.applyDDK(obj.C, obj.S, W);
        out.sigmaC = nan(size(obj.C)); out.sigmaS = nan(size(obj.S));
        out.history(end+1) = sprintf("DDK filter applied (%s, nmax=%d)", ...
            W.name, W.nmax);
    end

    % ---------------------------------------------------------- gfct, I/O
    function out = evalAt(obj, epoch)
        %EVALAT Evaluate gfct time-variable terms at an epoch.
        %   Outputs
        %     out        (1 x 1) shCoefficients   gfct variable terms evaluated at the epoch; static models pass through unchanged
        %
        %   Outputs
        %     out  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            epoch (1,1) double
        end
        model = struct('C', obj.C, 'S', obj.S, ...
            'variableTerms', obj.variableTerms);
        out = obj;
        [out.C, out.S] = shLowLevel.shEvalGFCT(model, epoch);
        out.epoch = epoch;
        out.history(end+1) = sprintf("gfct evaluated at %.4f", epoch);
    end

    function write(obj, filename, opts)
        %WRITE Export to an ICGEM .gfc file (round-trip tested).
        %
        %
        %   Inputs
        %     filename  char/string  path of the file to read/write (gzipped .gz accepted where documented)
        %   Options
        %     Comment ("")  free-text comment written into the gfc header
        %     Sidecar (true)  also write "<filename>.provenance.json".
        %             Set false when exporting into a folder that will be
        %             read back as a series: the sidecar is metadata, and
        %             a loose read pattern would otherwise pick it up
        arguments
            obj
            filename {mustBeTextScalar}
            opts.Comment (1,1) string = ""
            opts.Sidecar (1,1) logical = true
        end
        cmt = strjoin(obj.history, " | ");
        if strlength(opts.Comment) > 0, cmt = opts.Comment + " | " + cmt; end
        nm = obj.name; if strlength(nm) == 0, nm = "shAnalysis_export"; end
        shLowLevel.writeGFC(filename, obj.C, obj.S, obj.GM, obj.R, ...
            SigmaC = obj.sigmaC, SigmaS = obj.sigmaS, ...
            ModelName = nm, TideSystem = obj.tideSystem, Comment = cmt, ...
            Sidecar = opts.Sidecar);
    end

    % ---------------------------------------------------- tvANS interface
    function x = vec(obj, idx)
        %VEC Pack into a coefficient vector in shLowLevel.shIndex ordering.
        %   Outputs
        %     x          (P x 1) double   coefficients in idx ordering
        %
        %   Outputs
        %     x    (P x 1) double  stacked coefficients in shLowLevel.shIndex order
        %     sig  (P x 1) double  matching sigmas ([] when absent)
        arguments
            obj
            idx (1,1) struct
        end
        if idx.Lmax > obj.nmax
            error('shCoefficients:badIndex', ...
                'idx.Lmax=%d exceeds nmax=%d.', idx.Lmax, obj.nmax);
        end
        x = shLowLevel.vecFromCS(obj.C, obj.S, idx);
    end

    function disp(obj)
        %DISP Compact display: size, epoch, GM/R, processing history.
        %   disp(G) prints one header line with the model name, one
        %   contract line (nmax, epoch, product type, tide system, GM,
        %   R), whether formal sigmas are attached and how many gfct
        %   variable terms are carried, followed by the processing
        %   history. Called automatically when an shCoefficients object
        %   is shown without a semicolon.
        %
        %   Inputs
        %     obj  (1,1) shCoefficients  the coefficient set to display
        %   Outputs
        %     none - the summary is written to the command window
        %
        %   Example
        %     g = shCoefficients.read("ITSG-Grace2018_n60_2008-04.gfc");
        %     disp(g)                            % degree, epoch, history
        %
        %   Developed by Matthias Weigelt with the help of Claude (Fable 5).
        fprintf('  shCoefficients: %s\n', obj.name);
        fprintf('    nmax=%d | epoch=%.4f | product=%s | tide=%s | GM=%.6e | R=%.1f\n', ...
            obj.nmax, obj.epoch, obj.productType, obj.tideSystem, obj.GM, obj.R);
        fprintf('    sigmas: %s | variableTerms: %d\n', ...
            string(~isempty(obj.sigmaC)), numel(obj.variableTerms));
        fprintf('    history:\n');
        fprintf('      %s\n', obj.history);
    end
end

methods (Static)
    function obj = read(filename, opts)
        %READ Read an ICGEM .gfc/.gfct/.gfc.gz file into an shCoefficients.
        %   G = shCoefficients.read(FILE) parses product type and mid-epoch
        %   from standard GRACE L2 filenames (GSM/GAA/GAB/GAC/GAD-2_...);
        %   override with Epoch (NaN)=..., ProductType ("")=... when needed.
        %
        %   Inputs
        %     filename  char/string  path of the file to read/write (gzipped .gz accepted where documented)
        %   Outputs
        %     obj        (1 x 1) shCoefficients   parsed gfc/gfct(.gz) with GM, R, sigmas, variable terms, Epoch= as given
        %
        %   Outputs
        %     obj  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            filename {mustBeTextScalar}
            opts.Epoch (1,1) double = NaN
            opts.ProductType (1,1) string = ""
        end
        m = shLowLevel.shReadGFC(char(filename));
        meta = shLowLevel.parseGraceFilename(filename);
        ep = opts.Epoch;
        if isnan(ep), ep = meta.epoch; end
        pt = opts.ProductType;
        if strlength(pt) == 0
            pt = meta.productType;
            if pt == "unknown" && isempty(m.variableTerms) ...
                    && ~isfinite(ep)
                pt = "static";
            end
        end
        ts = "unknown"; nm = "";
        if isfield(m.header, 'tide_system'), ts = string(m.header.tide_system); end
        if isfield(m.header, 'modelname'),   nm = string(m.header.modelname); end
        [~, base, ext] = fileparts(char(filename));
        if strlength(nm) == 0
            % headers without modelname (COST-G GSM among them) used to
            % leave the name empty all the way into shSeries.names; the
            % filename stem is always available and always identifying
            nm = string(base);
        end
        obj = shCoefficients(m.C, m.S, SigmaC = m.sigmaC, SigmaS = m.sigmaS, ...
            GM = m.GM, R = m.R, Epoch = ep, ProductType = pt, ...
            TideSystem = ts, Name = nm, Header = m.header, ...
            VariableTerms = m.variableTerms, ...
            History = sprintf("read %s%s (nmax=%d)", base, ext, m.nmax));
    end

    function obj = fromVec(x, idx, like)
        %FROMVEC Unpack a shLowLevel.shIndex-ordered vector into an shCoefficients.
        %   G = shCoefficients.fromVec(X, IDX, LIKE) copies GM/R and other
        %   metadata from the template LIKE (an shCoefficients).
        %   Outputs
        %     obj        (1 x 1) shCoefficients   rebuilt from x with sizes/GM/R/epoch of the template
        %
        %   Outputs
        %     obj  (1,1) shCoefficients  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            x (:,1) double
            idx (1,1) struct
            like (1,1) shCoefficients
        end
        [cnm, snm] = shLowLevel.csFromVec(x, idx);
        obj = shCoefficients(cnm, snm, GM = like.GM, R = like.R, ...
            Epoch = like.epoch, ProductType = like.productType, ...
            TideSystem = like.tideSystem, Name = like.name + "_fromVec", ...
            History = [like.history; "rebuilt from vector"]);
    end

    function [obj, info] = analysis(grid, latVec, lonVec, nmax, opts)
        %ANALYSIS Spherical harmonic analysis: gridded/scattered data ->
        %   Stokes coefficients (the inverse of synthesis).
        %
        %   [OBJ, INFO] = shCoefficients.analysis(GRID, LAT, LON, NMAX,
        %       quantity ("geoid")=..., Method ("auto")=..., Weights ("none")=..., Kaula (0)=..., ...)
        %   estimates an shCoefficients object from a ring grid
        %   (GRID nlat x nlon, uniform full-circle LON: exact fast
        %   per-order solver) or from scattered points (equal-length
        %   GRID/LAT/LON vectors: full least squares). See
        %   shLowLevel.shAnalysisGrid for the estimation options, the Kaula
        %   regularization needed for under-determined sampling, and the
        %   exactness guarantees. GM (3.986004415e14)/R (6378136.3)/Epoch (NaN)/Name ("analysis") flow into the object.
        %
        %   Inputs   grid (nlat,nlon) or (Np,1) double, lat/lon [deg]
        %     lonVec  (P,1) double  observation longitudes [deg]
        %     latVec  (P,1) double  observation latitudes [deg, geocentric]
        %            nmax (1,1) double
        %   Outputs  obj  shCoefficients;  info: see shLowLevel.shAnalysisGrid
        %   Outputs
        %     obj        (1 x 1) shCoefficients   Stokes coefficients estimated from the grid (exact on ring grids; Kaula for scattered points)
        %
        %   Options
        %     kn ([])  load Love numbers, degrees 0..nmax (user-supplied; e.g. shLowLevel.fetchLoveNumbers)
        %     hn ([])  vertical-deformation Love numbers, degrees 0..nmax (user-supplied)
        %     rho_ave (5517)  mean Earth density [kg/m^3] in the EWH kernel
        %     rho_water (1000)  water density [kg/m^3] in the EWH kernel
        %     LatType ("geocentric")  "geocentric" | "geodetic": how the
        %         latitude inputs are interpreted; geodetic values are
        %         converted with Flattening before evaluation
        %     Flattening (1/298.257223563)  flattening of the reference
        %         ellipsoid for the geodetic/geocentric conversion (WGS84;
        %         overridable, never silently assumed)
        arguments
            grid double
            latVec double
            lonVec double
            nmax (1,1) double
            opts.quantity (1,1) string = "geoid"
            opts.Method (1,1) string = "auto"
            opts.Weights (1,1) string = "none"
            opts.Kaula (1,1) double = 0
            opts.kn double = []
            opts.hn double = []
            opts.rho_ave (1,1) double = 5517
            opts.rho_water (1,1) double = 1000
            opts.GM (1,1) double = 3.986004415e14
            opts.R (1,1) double = 6378136.3
            opts.Epoch (1,1) double = NaN
            opts.Name (1,1) string = "analysis"
            opts.LatType (1,1) string ...
                {mustBeMember(opts.LatType, ["geocentric","geodetic"])} = "geocentric"
            opts.Flattening (1,1) double = 1/298.257223563
        end
        if opts.LatType == "geodetic"
            latVec = shLowLevel.geodetic2geocentric(latVec, ...
                Flattening = opts.Flattening);
        end
        [C, S, info] = shLowLevel.shAnalysisGrid(grid, latVec, lonVec, nmax, ...
            Method = opts.Method, quantity = char(opts.quantity), ...
            GM = opts.GM, R = opts.R, kn = opts.kn, hn = opts.hn, ...
            rho_ave = opts.rho_ave, rho_water = opts.rho_water, ...
            Weights = opts.Weights, Kaula = opts.Kaula);
        obj = shCoefficients(C, S, GM = opts.GM, R = opts.R, ...
            Epoch = opts.Epoch, Name = opts.Name, ProductType = "analysis", ...
            History = sprintf("SH analysis (%s, %d points, residRMS %.3g)", ...
                info.method, info.nPoints, info.residRMS));
    end
end

methods (Static, Access = private)
    function [a, b] = checkPair(a, b, opname)
        if ~isa(a, 'shCoefficients') || ~isa(b, 'shCoefficients')
            error('shCoefficients:badInput', ...
                '%s requires two shCoefficients objects.', opname);
        end
        if a.nmax ~= b.nmax
            error('shCoefficients:sizeMismatch', ...
                'nmax differs (%d vs %d) - use truncate() first.', ...
                a.nmax, b.nmax);
        end
        if abs(a.GM - b.GM) > 1e-9 * a.GM || abs(a.R - b.R) > 1e-9 * a.R
            error('shCoefficients:constantsMismatch', ...
                ['GM/R differ between operands - convert first with ' ...
                 'toReference(GM=..., R=...) (exact spectral rescaling).']);
        end
    end

    function checkEpochs(a, b, tol, opname)
        if isfinite(a.epoch) && isfinite(b.epoch) && abs(a.epoch - b.epoch) > tol
            error('shCoefficients:epochMismatch', ...
                '%s: epochs %.4f and %.4f differ by more than %.3f yr.', ...
                opname, a.epoch, b.epoch, tol);
        end
    end

    function t = combinedType(ta, tb, sym)
        gax = ["GAA","GAB","GAC","GAD"];
        if ta == "GSM" && any(tb == gax)
            t = "GSM" + sym + tb;
        elseif tb == "GSM" && any(ta == gax)
            t = "GSM" + sym + ta;
        elseif ta == tb
            t = ta;
        else
            t = "combined";
        end
    end
end

methods (Access = private)
    function out = combine(a, b, sgn, newType)
        out = a;
        out.C = a.C + sgn * b.C;
        out.S = a.S + sgn * b.S;
        if ~isempty(a.sigmaC) && ~isempty(b.sigmaC)
            out.sigmaC = sqrt(a.sigmaC.^2 + b.sigmaC.^2);
        else
            out.sigmaC = [];
        end
        if ~isempty(a.sigmaS) && ~isempty(b.sigmaS)
            out.sigmaS = sqrt(a.sigmaS.^2 + b.sigmaS.^2);
        else
            out.sigmaS = [];
        end
        out.productType = newType;
        out.variableTerms = [];
        if sgn > 0
            out.history(end+1) = sprintf("plus: %s + %s -> %s", ...
                a.productType, b.productType, newType);
        else
            out.history(end+1) = "minus";
        end
    end
end
end
