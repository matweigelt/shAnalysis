classdef shSeries
%SHSERIES Time series of spherical harmonic coefficient sets.
%
%   Wraps a stack of shCoefficients over epochs and provides everything
%   with a time dimension: the mean field, the climatology (bias / trend /
%   annual / semi-annual), GAX background restoration, the classic
%   destripe+Gaussian chain applied per epoch, and the tvANS
%   time-variable anisotropic Wiener filter with exact basin
%   deconvolution. Immutable value class with provenance in .history.
%
%   Construction
%     ts = shSeries.read("GSM-2_*.gfc");        % pattern, sorted by epoch
%     ts = shSeries.read(fileList);              % string array of paths
%     ts = shSeries(objArray);                   % from shCoefficients array
%
%   Typical use
%     ts   = ts.restore(shSeries.read("GAD-2_*.gfc"));   % add background
%     m    = ts.mean;                                    % mean field
%     dts  = ts - m;                                     % anomalies
%     [clim, resid] = dts.climatology(Robust=true);      % + shClimatology
%     tsA  = dts.destripe(minOrder=6).gaussian(300);     % classic chain
%     [tsF, op] = dts.filter("tvANS", Constraints=ocean);% optimal chain
%     avg  = tsF.basinAverage(B, Deconvolve=true, Op=op);
%
%   Properties (read-only)
%     Cs, Ss        (nmax+1)x(nmax+1)xT double   coefficient stacks
%     sigmaCs/Ss    same or []                    formal error stacks
%     epochs        (T,1) double  decimal years
%     GM, R         (1,1) double
%     nmax, nEpochs (dependent)
%     productType   (1,1) string
%     names         (T,1) string
%     history       string column
%
%   See also shCoefficients, shClimatology.
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

properties (SetAccess = private)
    Cs double
    Ss double
    sigmaCs double = []
    sigmaSs double = []
    epochs (:,1) double
    GM (1,1) double = 3.986004415e14
    R (1,1) double = 6378136.3
    productType (1,1) string = "unknown"
    names (:,1) string = string.empty(0,1)
    history string = string.empty(0,1)
end

properties (Dependent)
    nmax
    nEpochs
end

methods
    function obj = set.history(obj, h)
        % normalize to a string column: scalar-seeded (end+1) growth
        % otherwise produces rows and breaks vertcat chains (v2.2 fix)
        h = string(h);                       % chars -> scalar string first
        obj.history = h(:);
    end

    function obj = shSeries(in, opts)
        %SHSERIES Construct from an shCoefficients array or raw stacks.
        %   ts = shSeries(objArray)
        %   ts = shSeries(Cs, Epochs ([])=..., Ss ([])=...)  (raw-stack form via opts)
        %
        %   Inputs
        %     in  (nmax+1 x nmax+1 x T) double or (1,T) shCoefficients  coefficient stacks or object array
        %   Outputs
        %     obj        (1 x 1) shSeries   epoch-sorted monthly stack with sigmas and history
        %
        %   Options
        %     GM (3.986004415e14)  gravitational constant times mass
        %         [m^3/s^2] of the series (overridable default)
        %     R (6378136.3)  reference radius [m] of the series
        %     ProductType ("unknown")  provider product code shared by all
        %         epochs, e.g. "GSM"
        %     Names (string.empty(0,1))  display labels, one per solution/series
        %     SigmaCs ([])  formal 1-sigma stack for Cs, same size
        %     SigmaSs ([])  formal 1-sigma stack for Ss, same size
        %     History (string.empty(0,1))  initial processing history; every
        %         operation appends one line
        %
        %   Outputs
        %     obj  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            in
            opts.Ss double = []
            opts.Epochs (:,1) double = []
            opts.GM (1,1) double = 3.986004415e14
            opts.R (1,1) double = 6378136.3
            opts.ProductType (1,1) string = "unknown"
            opts.Names (:,1) string = string.empty(0,1)
            opts.SigmaCs double = []
            opts.SigmaSs double = []
            opts.History string = string.empty(0,1)
        end
        if isa(in, 'shCoefficients')
            arr = in(:);
            T = numel(arr);
            nm = arr(1).nmax;
            for k = 2:T
                if arr(k).nmax ~= nm
                    error('shSeries:sizeMismatch', ...
                        'All epochs must share nmax (epoch %d has %d, expected %d) - truncate first.', ...
                        k, arr(k).nmax, nm);
                end
                if abs(arr(k).GM - arr(1).GM) > 1e-9*arr(1).GM || ...
                   abs(arr(k).R - arr(1).R) > 1e-9*arr(1).R
                    error('shSeries:constantsMismatch', ...
                        'GM/R differ within the series (epoch %d).', k);
                end
            end
            obj.Cs = zeros(nm+1, nm+1, T);
            obj.Ss = zeros(nm+1, nm+1, T);
            hasSig = all(arrayfun(@(g) ~isempty(g.sigmaC), arr));
            if hasSig
                obj.sigmaCs = zeros(nm+1, nm+1, T);
                obj.sigmaSs = zeros(nm+1, nm+1, T);
            end
            ep = zeros(T,1); nmv = strings(T,1);
            for k = 1:T
                obj.Cs(:,:,k) = arr(k).C;
                obj.Ss(:,:,k) = arr(k).S;
                if hasSig
                    obj.sigmaCs(:,:,k) = arr(k).sigmaC;
                    obj.sigmaSs(:,:,k) = arr(k).sigmaS;
                end
                ep(k) = arr(k).epoch;
                nmv(k) = arr(k).name;
            end
            [ep, ord] = sortWithNaN(ep);
            obj.Cs = obj.Cs(:,:,ord); obj.Ss = obj.Ss(:,:,ord);
            if hasSig
                obj.sigmaCs = obj.sigmaCs(:,:,ord);
                obj.sigmaSs = obj.sigmaSs(:,:,ord);
            end
            obj.epochs = ep;
            obj.names = nmv(ord);
            obj.GM = arr(1).GM; obj.R = arr(1).R;
            obj.productType = arr(1).productType;
            obj.history = sprintf("series from %d shCoefficients (nmax=%d)", T, nm);
        else
            % raw-stack form
            obj.Cs = in;
            obj.Ss = opts.Ss;
            if isempty(obj.Ss) || ~isequal(size(obj.Cs), size(obj.Ss)) ...
                    || size(in,1) ~= size(in,2)
                error('shSeries:badInput', ...
                    'Raw-stack form requires square Cs and matching Ss=... stacks.');
            end
            T = size(in, 3);
            if isempty(opts.Epochs)
                error('shSeries:badInput', 'Raw-stack form requires Epochs=...');
            end
            if numel(opts.Epochs) ~= T
                error('shSeries:badInput', 'numel(Epochs) must equal size(Cs,3).');
            end
            obj.epochs = opts.Epochs;
            obj.GM = opts.GM; obj.R = opts.R;
            obj.productType = opts.ProductType;
            obj.names = opts.Names;
            obj.sigmaCs = opts.SigmaCs; obj.sigmaSs = opts.SigmaSs;
            if isempty(opts.History)
                obj.history = sprintf("series from raw stacks (T=%d)", T);
            else
                obj.history = opts.History;
            end
        end
    end

    function n = get.nmax(obj), n = size(obj.Cs, 1) - 1; end
    function T = get.nEpochs(obj), T = size(obj.Cs, 3); end

    function g = at(obj, k)
        %AT Extract epoch k as an shCoefficients.
        %   Outputs
        %     g          (1 x 1) shCoefficients   month k with epoch and sigmas
        %
        %   Outputs
        %     g  (1,1) shCoefficients  the k-th epoch as a standalone field
        arguments
            obj
            k (1,1) double {mustBeInteger, mustBePositive}
        end
        if k > obj.nEpochs
            error('shSeries:badIndex', 'k=%d exceeds nEpochs=%d.', k, obj.nEpochs);
        end
        sc = []; ss = [];
        if ~isempty(obj.sigmaCs), sc = obj.sigmaCs(:,:,k); ss = obj.sigmaSs(:,:,k); end
        nm = "";
        if numel(obj.names) >= k, nm = obj.names(k); end
        g = shCoefficients(obj.Cs(:,:,k), obj.Ss(:,:,k), SigmaC = sc, ...
            SigmaS = ss, GM = obj.GM, R = obj.R, Epoch = obj.epochs(k), ...
            ProductType = obj.productType, Name = nm, ...
            History = [obj.history; sprintf("extracted epoch %d", k)]);
    end

    % --------------------------------------------------- series statistics
    function g = mean(obj, opts)
        %MEAN Mean field of the series.
        %   G = ts.mean (Omitnan=true default) returns an shCoefficients;
        %   its sigmas, when present, are the standard error of the mean.
        %   Outputs
        %     g          (1 x 1) shCoefficients   temporal mean field (omits NaN months)
        %
        %   Outputs
        %     out  (1,1) shCoefficients  epoch-mean field (sigmas RSS/T when present)
        arguments
            obj
            opts.Omitnan (1,1) logical = true
        end
        if opts.Omitnan, flag = 'omitnan'; else, flag = 'includenan'; end
        Cm = mean(obj.Cs, 3, flag);
        Sm = mean(obj.Ss, 3, flag);
        T = obj.nEpochs;
        sc = std(obj.Cs, 0, 3, flag) / sqrt(T);
        ss = std(obj.Ss, 0, 3, flag) / sqrt(T);
        g = shCoefficients(Cm, Sm, SigmaC = sc, SigmaS = ss, GM = obj.GM, ...
            R = obj.R, Epoch = mean(obj.epochs, 'omitnan'), ...
            ProductType = obj.productType, Name = "mean(" + obj.productType + ")", ...
            History = [obj.history; sprintf("mean of %d epochs", T)]);
    end

    function [clim, resid] = climatology(obj, opts)
        %CLIMATOLOGY Bias/trend/annual/semi-annual (+extra periods) fit.
        %   [CLIM, RESID] = ts.climatology(Robust=false, T0 (NaN)=mean(epochs),
        %       Periods (double.empty(1,0))=[], Weights ([])=[])
        %   CLIM is an shClimatology; RESID the residual shSeries.
        %   Periods: extra periodic terms [years] - e.g. the GRACE tidal
        %   alias periods [161/365.25, 3.66, 7.48] (S2, K2, K1); without
        %   them, trend and semi-annual estimates absorb tidal aliasing.
        %   Weights: T x 1 per-epoch weights (e.g. 1./info.vce from the
        %   tvANS filter). CLIM carries per-coefficient 1-sigma
        %   significance (OLS formula, Monte-Carlo validated); component
        %   accessors return shCoefficients WITH sigmas, enabling
        %   significance-masked trend maps. Requires >= 6+2K epochs and
        %   NaN-free stacks (use dropNaN/select across the 2017-2018 gap).
        %   Outputs
        %     clim       (1 x 1) shClimatology   fitted bias/trend/annual/semiannual (+Periods=) with coefficient sigmas
        %     resid      (1 x 1) shSeries   residual series about the fit
        %
        %   Options
        %     ARCorrect (false)  correct the parameter sigmas for AR(1)
        %         residual autocorrelation (Kendall-corrected r1); the
        %         uncorrected sigmas are optimistic for monthly series
        arguments
            obj
            opts.Robust (1,1) logical = false
            opts.T0 (1,1) double = NaN
            opts.Periods (1,:) double = double.empty(1,0)
            opts.Weights (:,1) double = []
            opts.ARCorrect (1,1) logical = false
        end
        obj.assertClean('climatology');
        if any(~isfinite(obj.epochs))
            error('shSeries:noEpoch', ...
                'climatology requires finite epochs for every entry.');
        end
        nPar = 6 + 2*numel(opts.Periods);
        if obj.nEpochs < nPar
            error('shSeries:tooFewEpochs', ...
                'climatology needs >= %d epochs for %d parameters (got %d).', ...
                nPar, nPar, obj.nEpochs);
        end
        t0 = opts.T0;
        if isnan(t0), t0 = mean(obj.epochs); end
        X = obj.flatten();                             % (2*Nc) x T
        [~, XresX, coef, ~, coefSigma] = shLowLevel.fitDeterministicModel( ...
            X, obj.epochs, Robust = opts.Robust, T0 = t0, ...
            Periods = opts.Periods, Weights = opts.Weights, ...
            ARCorrect = opts.ARCorrect);
        clim = shClimatology.fromCoef(coef, t0, obj, ...
            Periods = opts.Periods, CoefSigma = coefSigma);
        resid = obj;
        [resid.Cs, resid.Ss] = obj.unflatten(XresX);
        if opts.ARCorrect
            clim = clim.withNote( ...
                "coefficient sigmas AR(1)-corrected (sqrt((1+r1)/(1-r1)))");
        end
        if isempty(opts.Periods)
            resid.history(end+1) = "climatology removed (bias/trend/annual/semi-annual)";
        else
            resid.history(end+1) = sprintf( ...
                "climatology removed (+%d extra periods)", numel(opts.Periods));
        end
    end

    function out = trendBreaks(obj, opts)
        %TRENDBREAKS Piecewise-linear trends with break significance.
        %   OUT = ts.trendBreaks(Breaks=[2016.0], Periods (double.empty(1,0))=[], T0 (NaN)=,
        %       ARCorrect (false)=false) fits the standard climatology design
        %   PLUS hinge terms max(0, t - tb) per break epoch (continuous
        %   piecewise-linear trend) to every coefficient, and tests the
        %   break's significance with an F-test against the no-break
        %   model (p-values via betainc - base MATLAB, no toolbox).
        %   Python-validated: 5.1% null rejection at the 5% level;
        %   F = 23 for a 0.02/yr break in 0.05-sigma noise (T = 120).
        %
        %   Returns a struct (NOT an shClimatology - the climatology
        %   schema stays fixed): trend1C/S (pre-break rate), hingeC/S
        %   (n1 x n1 x K rate CHANGES per break; post-break rate =
        %   trend1 + cumsum of hinges), sigmas, F (n1 x n1 x K),
        %   pValue, breaks, t0.
        %   Outputs
        %     out        struct: trends per segment (shCoefficients), F/p per break, segment epochs
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            opts.Breaks (1,:) double
            opts.Periods (1,:) double = double.empty(1,0)
            opts.T0 (1,1) double = NaN
            opts.ARCorrect (1,1) logical = false
        end
        assert(~isempty(opts.Breaks), 'shSeries:noBreaks', ...
            'trendBreaks needs at least one break epoch.');
        t0 = opts.T0;
        if isnan(t0), t0 = mean(obj.epochs); end
        n1 = obj.nmax + 1; Nc = n1^2; T = obj.nEpochs;
        X = [reshape(obj.Cs, Nc, T); reshape(obj.Ss, Nc, T)];
        K = numel(opts.Breaks);
        nP = numel(opts.Periods);
        [~, Xr1, c1, ~, cs1] = shLowLevel.fitDeterministicModel(X, obj.epochs, ...
            T0 = t0, Periods = opts.Periods, ARCorrect = opts.ARCorrect);
        [~, XrB, cB, ~, csB] = shLowLevel.fitDeterministicModel(X, obj.epochs, ...
            T0 = t0, Periods = opts.Periods, Breaks = opts.Breaks, ...
            ARCorrect = opts.ARCorrect);
        ss0 = sum(Xr1.^2, 2); ss1 = sum(XrB.^2, 2);
        d2 = T - (6 + 2*nP + K);
        F = max((ss0 - ss1) / K, 0) ./ max(ss1 / max(d2, 1), realmin);
        % p = 1 - Fcdf(F; K, d2) via the regularized incomplete beta
        pv = betainc(d2 ./ (d2 + K * F), d2/2, K/2);
        rowH = 6 + 2*nP + (1:K);
        out = struct('breaks', opts.Breaks, 't0', t0, ...
            'trend1C', reshape(cB(2, 1:Nc), n1, n1), ...
            'trend1S', reshape(cB(2, Nc+1:end), n1, n1), ...
            'sigmaTrend1C', reshape(csB(2, 1:Nc), n1, n1), ...
            'sigmaTrend1S', reshape(csB(2, Nc+1:end), n1, n1), ...
            'trendNoBreakC', reshape(c1(2, 1:Nc), n1, n1), ...
            'sigmaTrendNoBreakC', reshape(cs1(2, 1:Nc), n1, n1), ...
            'hingeC', zeros(n1, n1, K), 'hingeS', zeros(n1, n1, K), ...
            'sigmaHingeC', zeros(n1, n1, K), 'sigmaHingeS', zeros(n1, n1, K), ...
            'F', zeros(n1, n1), 'pValue', zeros(n1, n1));
        for k = 1:K
            out.hingeC(:, :, k) = reshape(cB(rowH(k), 1:Nc), n1, n1);
            out.hingeS(:, :, k) = reshape(cB(rowH(k), Nc+1:end), n1, n1);
            out.sigmaHingeC(:, :, k) = reshape(csB(rowH(k), 1:Nc), n1, n1);
            out.sigmaHingeS(:, :, k) = reshape(csB(rowH(k), Nc+1:end), n1, n1);
        end
        out.F = reshape(F(1:Nc), n1, n1);
        out.pValue = reshape(pv(1:Nc), n1, n1);
        out.FS = reshape(F(Nc+1:end), n1, n1);
        out.pValueS = reshape(pv(Nc+1:end), n1, n1);
    end

    function out = fan(obj, rDegKm, rOrdKm)
        %FAN Han fan filter on every epoch (degree x order Gaussian).
        %   OUT = ts.fan(RDEG, RORD)  - see shLowLevel.shFanFilter.
        %
        %   Inputs
        %     rOrdKm  (1,1) double  fan filter order-direction radius [km]
        %     rDegKm  (1,1) double  fan filter degree-direction radius [km]
        %   Outputs
        %     out        (1 x 1) shSeries   fan-filtered per month
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            rDegKm (1,1) double {mustBePositive}
            rOrdKm (1,1) double {mustBePositive}
        end
        out = obj;
        for k = 1:obj.nEpochs
            [out.Cs(:,:,k), out.Ss(:,:,k)] = shLowLevel.shFanFilter( ...
                obj.Cs(:,:,k), obj.Ss(:,:,k), rDegKm, rOrdKm, R = obj.R);
        end
        if ~isempty(out.sigmaCs)
            Wn = shLowLevel.shGaussianWeights(obj.nmax, rDegKm, R = obj.R/1e3);
            Wm = shLowLevel.shGaussianWeights(obj.nmax, rOrdKm, R = obj.R/1e3);
            Wmat = Wn(:) .* Wm(:)';
            for k = 1:obj.nEpochs
                out.sigmaCs(:,:,k) = obj.sigmaCs(:,:,k) .* Wmat;
                out.sigmaSs(:,:,k) = obj.sigmaSs(:,:,k) .* Wmat;
            end
        end
        out.history(end+1) = sprintf( ...
            "fan filtered (deg %g km, ord %g km)", rDegKm, rOrdKm);
    end

    function out = applyTN14(obj, tn, opts)
        %APPLYTN14 Replace C20 (and flagged C30) per month from TN-14.
        %   OUT = ts.applyTN14("TN-14.txt") - series-level convenience
        %   around shCoefficients.applyTN14 (v2.4.1; the per-month epoch
        %   selects the TN entry). All options pass through.
        %   Outputs
        %     out        (1 x 1) shSeries   C20/C30 replaced epoch-matched
        %
        %   Options
        %     Tolerance (0.05)  maximum |epoch difference| [yr] accepted
        %         when matching a table entry to this epoch
        %     ReplaceC30 ("auto")  "auto" (replace when the table has a
        %         non-NaN C30 for that month) | "never" | "always"
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            tn
            opts.Tolerance (1,1) double = 0.05
            opts.ReplaceC30 (1,1) string = "auto"
        end
        arr = shCoefficients.empty(0, 1);
        for k = 1:obj.nEpochs
            g = obj.at(k);
            arr(k, 1) = g.applyTN14(tn, Tolerance = opts.Tolerance, ...
                ReplaceC30 = opts.ReplaceC30);
        end
        out = shSeries(arr);
    end

    function out = addDegree1(obj, tn, opts)
        %ADDDEGREE1 Insert TN-13 geocenter degree-1 per month.
        %   OUT = ts.addDegree1("TN-13.txt") - series-level convenience
        %   around shCoefficients.addDegree1 (v2.4.1).
        %   Outputs
        %     out        (1 x 1) shSeries   degree 1 completed epoch-matched
        %
        %   Options
        %     Tolerance (0.05)  maximum |epoch difference| [yr] accepted
        %         when matching a table entry to this epoch
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            tn
            opts.Tolerance (1,1) double = 0.05
        end
        arr = shCoefficients.empty(0, 1);
        for k = 1:obj.nEpochs
            g = obj.at(k);
            arr(k, 1) = g.addDegree1(tn, Tolerance = opts.Tolerance);
        end
        out = shSeries(arr);
    end

    function out = removeAlias(obj, opts)
        %REMOVEALIAS Remove a tidal alias harmonic from every epoch.
        %   TS2 = ts.removeAlias() fits, per coefficient, a harmonic at
        %   the S2 tidal alias period (161 days) TOGETHER WITH bias,
        %   trend, annual and semi-annual, and subtracts only that
        %   harmonic. This is the correction GravIS applies to its
        %   Level-2B products (https://gravis.gfz.de/corrections):
        %   errors in the ocean-tide background model alias into the
        %   monthly fields at this period, and leaving them in puts a
        %   spurious ~161-day signal into every derived series.
        %
        %   Fitting the alias jointly with the deterministic terms
        %   matters: fitted alone it absorbs part of the trend and the
        %   annual cycle, and subtracting that would damage both.
        %
        %   PHASE OFFSET. Across the GRACE / GRACE-FO boundary the nodal
        %   planes are not aligned, so one harmonic cannot describe both
        %   missions. Landerer et al. (2020) prescribe a fixed 100 degree
        %   offset for the later mission; PhaseOffset/SplitEpoch apply
        %   it. The offset is FIXED, so the model gains no free
        %   parameters - it only gains the ability to fit both spans.
        %   Ignoring it mis-estimates the amplitude badly (validated in
        %   tools/dev/validate_alias.py: 0.41 error on a unit-scale
        %   synthetic). On a series that lies entirely on one side of
        %   SplitEpoch the offset has no effect.
        %
        %   Inputs
        %     (none beyond the object)
        %
        %   Options
        %     Period (161/365.25)  alias period [years]. The default is
        %             the S2 alias; the K2 and K1 aliases (3.66 and 7.48
        %             years) are longer than most records and are better
        %             handled as climatology Periods
        %     PhaseOffset (100)  phase offset [degrees] applied to epochs
        %             at or after SplitEpoch
        %     SplitEpoch (2018.0)  mission boundary [decimal year]; the
        %             GRACE-FO record starts mid-2018
        %     T0 (NaN)  reference epoch of the fit (NaN: the mean epoch)
        %
        %   Outputs
        %     out        (1,1) shSeries  copy with the alias harmonic
        %                removed from every epoch; the operation is
        %                appended to the history
        %
        %   Example
        %     ts = ts.removeAlias();                   % S2, 100 deg offset
        %     ts = ts.removeAlias(SplitEpoch = 2018.4);
        %
        %   See also shSeries.climatology, shLowLevel.standardChain.
        %
        %   Developed by Matthias Weigelt with the help of Claude
        %   (Opus 5), 2026-08-11 (v3.4.0).
        arguments
            obj
            opts.Period (1,1) double {mustBePositive} = 161/365.25
            opts.PhaseOffset (1,1) double = 100
            opts.SplitEpoch (1,1) double = 2018.0
            opts.T0 (1,1) double = NaN
        end
        obj.assertClean('removeAlias');
        t = obj.epochs(:);
        if any(~isfinite(t))
            error('shSeries:noEpoch', ...
                'removeAlias requires finite epochs for every entry.');
        end
        t0 = opts.T0;
        if ~isfinite(t0), t0 = mean(t); end
        w = 2 * pi;
        a = w * t / opts.Period + deg2rad(opts.PhaseOffset) * ...
            (t >= opts.SplitEpoch);
        A = [ones(numel(t), 1), t - t0, cos(w*t), sin(w*t), ...
             cos(2*w*t), sin(2*w*t), cos(a), sin(a)];
        % A short record does NOT show up as rank deficiency or bad
        % conditioning - a one-year span gives a full-rank design with a
        % condition number of 38, and a fit that means nothing. The real
        % requirements are degrees of freedom and enough alias cycles to
        % separate the harmonic from the seasonal terms.
        nPar = size(A, 2);
        span = max(t) - min(t);
        if numel(t) < 2 * nPar || span < 2 * opts.Period || rank(A) < nPar
            error('shSeries:removeAlias:tooShort', ...
                ['%d epochs spanning %.2f yr cannot support the ' ...
                 '%d-term fit (bias, trend, annual, semi-annual, ' ...
                 'alias). At least %d epochs and %.2f yr (two alias ' ...
                 'periods) are required, and in practice several years ' ...
                 'are needed before the estimate is stable.'], ...
                numel(t), span, nPar, 2 * nPar, 2 * opts.Period);
        end
        n1 = obj.nmax + 1;
        Y = [reshape(obj.Cs, [], obj.nEpochs); ...
             reshape(obj.Ss, [], obj.nEpochs)].';      % T x 2*n1^2
        coef = A \ Y;
        alias = A(:, 7:8) * coef(7:8, :);              % ONLY the harmonic
        Z = (Y - alias).';
        out = obj;
        out.Cs = reshape(Z(1:n1^2, :), n1, n1, obj.nEpochs);
        out.Ss = reshape(Z(n1^2+1:end, :), n1, n1, obj.nEpochs);
        out.history(end+1) = sprintf( ...
            "alias removed (%.1f d, %.0f deg offset from %.2f)", ...
            opts.Period * 365.25, opts.PhaseOffset, opts.SplitEpoch);
    end

    function out = removeGIA(obj, gia, opts)
        %REMOVEGIA Subtract a glacial isostatic adjustment model.
        %   OUT = ts.removeGIA(GIA, T0 (NaN)=mean(epochs)) subtracts the linear
        %   GIA signal  giaRate * (t - T0)  from every epoch. GIA is an
        %   shCoefficients object holding RATE coefficients [1/yr] - read
        %   your model (ICE-6G_D, Caron et al., ...) from its gfc file
        %   with shCoefficients.read; models are ALWAYS user-supplied,
        %   never embedded. Trends computed after this are GIA-corrected;
        %   the model is treated as exact (its uncertainty, often the
        %   dominant trend error over Laurentia/Fennoscandia/Antarctica,
        %   is NOT propagated - documented limitation, see the guide).
        %   The GIA model is truncated to the series nmax; a model with
        %   smaller nmax is zero-padded with a note in history.
        %   Outputs
        %     out        (1 x 1) shSeries   GIA trend removed about the series mean epoch
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            gia (1,1) shCoefficients
            opts.T0 (1,1) double = NaN
        end
        t0 = opts.T0;
        if isnan(t0), t0 = mean(obj.epochs); end
        n1 = obj.nmax + 1;
        gC = zeros(n1); gS = zeros(n1);
        nc = min(n1, size(gia.C, 1));
        gC(1:nc, 1:nc) = gia.C(1:nc, 1:nc);
        gS(1:nc, 1:nc) = gia.S(1:nc, 1:nc);
        out = obj;
        for k = 1:obj.nEpochs
            dt = obj.epochs(k) - t0;
            out.Cs(:, :, k) = obj.Cs(:, :, k) - gC * dt;
            out.Ss(:, :, k) = obj.Ss(:, :, k) - gS * dt;
        end
        out.history(end+1) = sprintf( ...
            "GIA removed (model nmax=%d, T0=%.4f, rate applied as exact)", ...
            min(size(gia.C, 1) - 1, obj.nmax), t0);
        if nc < n1
            out.history(end+1) = sprintf( ...
                "GIA model zero-padded above degree %d", nc - 1);
        end
    end

    function out = applyDDK(obj, W)
        %APPLYDDK Apply a DDK decorrelation filter to every epoch.
        %   OUT = ts.applyDDK(W), W from shLowLevel.readDDK. Sigma stacks are
        %   invalidated (NaN) - the filter correlates coefficients.
        %   Outputs
        %     out        (1 x 1) shSeries   DDK-decorrelated per month
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            W
        end
        W = shLowLevel.readDDK(W);
        out = obj;
        for k = 1:obj.nEpochs
            [out.Cs(:, :, k), out.Ss(:, :, k)] = ...
                shLowLevel.applyDDK(obj.Cs(:, :, k), obj.Ss(:, :, k), W);
        end
        out.sigmaCs = []; out.sigmaSs = [];
        out.history(end+1) = sprintf("DDK filter applied (%s)", W.name);
    end

    % --------------------------------------------------------- GAX restore
    function [rep, h] = compare(obj, others, varargin)
        %COMPARE Temporal comparison against other series (v2.6.0).
        %   REP = ts.compare(TS2) or ts.compare({TS2, TS3, ...})
        %   delegates to shLowLevel.compareSeries with OBJ as the reference;
        %   all options pass through (Basin=, Names=, Plot=, ...).
        %
        %   Inputs
        %     others  shSeries or cell of shSeries  series to compare; OBJ is the reference
        %   Outputs
        %     rep        (1 x 1) struct   shLowLevel.compareSeries report
        %     h          (1 x 1) graphics handle   figure (Plot = true)
        %
        %   Outputs
        %     rep  (1,1) struct  shLowLevel.compareSeries report
        %     h    (1,1) graphics handle  figure (Plot = true only)
        if ~iscell(others), others = {others}; end
        [rep, h] = shLowLevel.compareSeries([{obj}, others], varargin{:});
    end

    function out = restore(obj, gax, opts)
        %RESTORE Add a background (GAA/GAB/GAC/GAD) series, epoch-matched.
        %   OUT = ts.restore(GAXSERIES, Tolerance=0.05, AllowMissing=false)
        %   Each epoch of OBJ is matched to the nearest GAX epoch within
        %   Tolerance [yr]. Missing matches error (or, with
        %   AllowMissing=true, leave those epochs unchanged with a note).
        %   Outputs
        %     out        (1 x 1) shSeries   background model added back epoch-matched (e.g. GSM + GAD)
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            gax (1,1) shSeries
            opts.Tolerance (1,1) double = 0.05
            opts.AllowMissing (1,1) logical = false
        end
        if obj.nmax ~= gax.nmax
            error('shSeries:sizeMismatch', ...
                'nmax differs (%d vs %d) - truncate first.', obj.nmax, gax.nmax);
        end
        out = obj;
        nMissing = 0;
        for k = 1:obj.nEpochs
            [d, j] = min(abs(gax.epochs - obj.epochs(k)));
            if isfinite(d) && d <= opts.Tolerance
                out.Cs(:,:,k) = obj.Cs(:,:,k) + gax.Cs(:,:,j);
                out.Ss(:,:,k) = obj.Ss(:,:,k) + gax.Ss(:,:,j);
            elseif opts.AllowMissing
                nMissing = nMissing + 1;
            else
                error('shSeries:epochMismatch', ...
                    'No %s epoch within %.3f yr of %.4f (epoch %d).', ...
                    gax.productType, opts.Tolerance, obj.epochs(k), k);
            end
        end
        out.productType = obj.productType + "+" + gax.productType;
        out.history(end+1) = sprintf("restored %s (%d/%d epochs matched)", ...
            gax.productType, obj.nEpochs - nMissing, obj.nEpochs);
    end

    % ----------------------------------------------------------- arithmetic
    function out = minus(a, b)
        %MINUS Subtract a field (shCoefficients) or an epoch-matched series.
        %
        %   Inputs
        %     b  (1,1) shSeries or shCoefficients  subtrahend (epoch-wise for series)
        %   Outputs
        %     out        (1 x 1) shSeries   per-month difference (series or single field)
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        if isa(b, 'shCoefficients')
            if a.nmax ~= b.nmax
                error('shSeries:sizeMismatch', ...
                    'nmax differs (%d vs %d).', a.nmax, b.nmax);
            end
            out = a;
            out.Cs = a.Cs - b.C;
            out.Ss = a.Ss - b.S;
            out.history(end+1) = "subtracted field: " + b.name;
        elseif isa(b, 'shSeries')
            if a.nmax ~= b.nmax
                error('shSeries:sizeMismatch', ...
                    'nmax differs (%d vs %d).', a.nmax, b.nmax);
            end
            if a.nEpochs ~= b.nEpochs || any(abs(a.epochs - b.epochs) > 0.05)
                error('shSeries:epochMismatch', ...
                    'Series minus requires matching epoch sets.');
            end
            out = a;
            out.Cs = a.Cs - b.Cs; out.Ss = a.Ss - b.Ss;
            out.history(end+1) = "subtracted series (" + b.productType + ")";
        else
            error('shSeries:badInput', 'Unsupported operand for minus.');
        end
    end

    % ------------------------------------------------------- classic chain
    function out = destripe(obj, opts)
        %DESTRIPE Swenson & Wahr / P3M6, applied per epoch.
        %   Outputs
        %     out        (1 x 1) shSeries   destriped per month
        %
        %   Options
        %     minOrder (6)  lowest order to destripe; coefficients at
        %         m < minOrder pass through unchanged (below it the signal
        %         is real, not stripes)
        %     polyOrder (3)  order of the polynomial fitted and removed
        %         per order/parity sequence
        %     windowLength ([])  [] fits ONE polynomial over the whole
        %         sequence; an odd integer >= polyOrder+2 uses a centered
        %         moving window instead (windowLength=6 with polyOrder=3 is
        %         the common "P3M6" variant)
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            opts.minOrder (1,1) double = 6
            opts.polyOrder (1,1) double = 3
            opts.windowLength double = []
        end
        obj.assertClean('destripe');
        out = obj;
        for k = 1:obj.nEpochs
            [out.Cs(:,:,k), out.Ss(:,:,k)] = shLowLevel.shDestripe( ...
                obj.Cs(:,:,k), obj.Ss(:,:,k), 'minOrder', opts.minOrder, ...
                'polyOrder', opts.polyOrder, 'windowLength', opts.windowLength);
        end
        out.history(end+1) = sprintf("destriped all epochs (minOrder=%d, polyOrder=%d)", ...
            opts.minOrder, opts.polyOrder);
    end

    function out = gaussian(obj, radiusKm)
        %GAUSSIAN Jekeli smoothing per epoch (vectorized over the stack).
        %
        %   Inputs
        %     radiusKm  (1,1) double  Gaussian filter half-response radius [km]
        %   Outputs
        %     out        (1 x 1) shSeries   Gaussian-smoothed per month
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            radiusKm (1,1) double {mustBePositive}
        end
        obj.assertClean('gaussian');
        Wn = shLowLevel.shGaussianWeights(obj.nmax, radiusKm);
        out = obj;
        out.Cs = obj.Cs .* Wn(:);
        out.Ss = obj.Ss .* Wn(:);
        if ~isempty(obj.sigmaCs)
            out.sigmaCs = obj.sigmaCs .* Wn(:);
            out.sigmaSs = obj.sigmaSs .* Wn(:);
        end
        out.history(end+1) = sprintf("Gaussian smoothed all epochs (%g km)", radiusKm);
    end

    function out = truncate(obj, nmaxNew)
        %TRUNCATE Reduce the maximum degree of the whole series.
        %
        %   Inputs
        %     nmaxNew  (1,1) double  new maximum degree (must not exceed the current nmax)
        %   Outputs
        %     out        (1 x 1) shSeries   truncated per month
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            nmaxNew (1,1) double {mustBeInteger, mustBeNonnegative}
        end
        if nmaxNew > obj.nmax
            error('shSeries:badTruncation', ...
                'Requested nmax %d exceeds available %d.', nmaxNew, obj.nmax);
        end
        out = obj; k = nmaxNew + 1;
        out.Cs = obj.Cs(1:k, 1:k, :); out.Ss = obj.Ss(1:k, 1:k, :);
        if ~isempty(obj.sigmaCs)
            out.sigmaCs = obj.sigmaCs(1:k,1:k,:);
            out.sigmaSs = obj.sigmaSs(1:k,1:k,:);
        end
        out.history(end+1) = sprintf("truncated to nmax=%d", nmaxNew);
    end

    % --------------------------------------------------- gap handling (v2.1)
    function out = dropNaN(obj)
        %DROPNAN Remove epochs containing any non-finite coefficient.
        %   OUT = ts.dropNaN() drops every epoch whose C or S stack has a
        %   NaN/Inf anywhere - typically placeholder months across the
        %   GRACE / GRACE-FO gap (2017.5-2018.5) or corrupt reads. The
        %   remaining series is irregularly sampled, which climatology and
        %   the tvANS filter handle natively (least-squares fit on the
        %   actual epochs; per-month VCE).
        %
        %   Outputs  out  shSeries (subset; unchanged if nothing to drop)
        %   Outputs
        %     out        (1 x 1) shSeries   gap months removed
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        bad = false(1, obj.nEpochs);
        for k = 1:obj.nEpochs
            bad(k) = any(~isfinite(obj.Cs(:,:,k)), 'all') ...
                  || any(~isfinite(obj.Ss(:,:,k)), 'all');
        end
        out = obj.select(~bad);
        out.history(end) = sprintf("dropNaN removed %d of %d epochs", ...
            nnz(bad), obj.nEpochs);
    end

    function out = select(obj, which)
        %SELECT Subset the series by index, mask, or time range.
        %   OUT = ts.select(WHICH) with WHICH one of
        %     logical (1,T)   keep masked epochs
        %     indices         keep listed epochs (in the given order)
        %     [tMin tMax]     keep epochs with tMin <= t <= tMax
        %                     (two-element non-integer-valued row => range;
        %                      use explicit logical/index forms otherwise)
        %
        %   Inputs   which  logical | double
        %   Outputs  out    shSeries subset
        %   Outputs
        %     out        (1 x 1) shSeries   months selected by logical/index mask
        %
        %   Outputs
        %     out  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            which {mustBeVector}
        end
        T = obj.nEpochs;
        if islogical(which)
            assert(numel(which) == T, 'shSeries:badInput', ...
                'Logical mask must have %d elements.', T);
            keep = find(which(:)');
        elseif numel(which) == 2 && ~all(which == round(which) & which >= 1 & which <= T)
            keep = find(obj.epochs(:)' >= which(1) & obj.epochs(:)' <= which(2));
        else
            assert(all(which == round(which)) && all(which >= 1) && all(which <= T), ...
                'shSeries:badInput', 'Indices must be integers in 1..%d.', T);
            keep = which(:)';
        end
        out = obj;
        out.Cs = obj.Cs(:,:,keep);
        out.Ss = obj.Ss(:,:,keep);
        if ~isempty(obj.sigmaCs), out.sigmaCs = obj.sigmaCs(:,:,keep); end
        if ~isempty(obj.sigmaSs), out.sigmaSs = obj.sigmaSs(:,:,keep); end
        out.epochs = obj.epochs(keep);
        if ~isempty(obj.names) && numel(obj.names) == T
            out.names = obj.names(keep);
        end
        out.history(end+1) = sprintf("select kept %d of %d epochs", ...
            numel(keep), T);
    end

    % ---------------------------------------------------------- tvANS chain
    function [out, op, info] = filter(obj, method, opts)
        %FILTER Optimal time-variable filtering of the series.
        %   [TSF, OP, INFO] = ts.filter("tvANS", Constraints ([])=..., NoiseCov ([])=...,
        %   SignalMode ("isotropic")="isotropic"|"inhomogeneous", NIterSignal (3)=3,
        %   Robust (false)=false, MinDegree (2)=2)
        %   runs the tvANS chain (see shLowLevel.tvANSFilter): deterministic-model
        %   separation, empirical or supplied noise covariance, monthly VCE
        %   scaling, data-driven signal covariance, per-month anisotropic
        %   Wiener filter, optional hard constraints. OP is the stored
        %   linear operator - required for basinAverage(...,
        %   Deconvolve=true) and shLowLevel.resolutionMap.
        %   v2.1: Blocks ("auto")="auto"|"on"|"off" engages the block-diagonal
        %   fast path (isotropic + empirical N + no constraints; identical
        %   results, tractable to Lmax ~ 120). The filtered series carries
        %   per-coefficient posterior 1-sigma stacks in OUT.sigmaCs/Ss
        %   (from info.sigmaXfres; degrees below MinDegree are NaN).
        %   v3.1.1: the three noise-covariance tuning options of
        %   shLowLevel.tvANSFilter are forwarded, so the class method is the
        %   full single point of access: Shrinkage ([])=empirical noise-
        %   covariance shrinkage, VCEMinDegree ([])=first degree entering
        %   the variance-component estimation, VCEBands ([])=order-band
        %   edges for per-band monthly VCE factors (block path only).
        %   [] means "leave the shLowLevel.tvANSFilter default in place"
        %   (0.1, round(2/3*Lmax) and no banding respectively) - the
        %   defaults are NOT duplicated here, they have exactly one home.
        %
        %   Inputs
        %     method  (1,1) string  filter chain to run; currently "tvANS"
        %             (the simple filters are their own methods:
        %             gaussian / fan / destripe / applyDDK)
        %   Outputs
        %     out        (1 x 1) shSeries   tvANS-filtered series with sigmaCs/Ss posterior stacks (exact incl. constraints, v2.5)
        %     op         struct   stored linear operator (V/Ut/lam/s or blocks, model, detLeverage/detResVar) for deconvolution and resolution maps
        %     info       struct: sigmaXfres (P x T), sigmaNote, VCE diagnostics
        %
        arguments
            obj
            method (1,1) string {mustBeMember(method, "tvANS")} %#ok<INUSA>
            opts.NoiseCov double = []
            opts.Constraints double = []
            opts.SignalMode (1,1) string = "isotropic"
            opts.NIterSignal (1,1) double = 3
            opts.Robust (1,1) logical = false
            opts.MinDegree (1,1) double = 2
            opts.Blocks (1,1) string ...
                {mustBeMember(opts.Blocks, ["auto","on","off"])} = "auto"
            opts.Shrinkage double {mustBeScalarOrEmpty} = []
            opts.VCEMinDegree double {mustBeScalarOrEmpty} = []
            opts.VCEBands (1,:) double = []
        end
        obj.assertClean('filter');
        if any(~isfinite(obj.epochs))
            error('shSeries:noEpoch', ...
                'filter requires finite epochs for every entry.');
        end
        idx = shLowLevel.shIndex(obj.nmax, MinDegree = opts.MinDegree);
        X = zeros(idx.P, obj.nEpochs);
        for k = 1:obj.nEpochs
            X(:,k) = shLowLevel.vecFromCS(obj.Cs(:,:,k), obj.Ss(:,:,k), idx);
        end
        % forward the tuning options only when the caller set them, so the
        % single home of their defaults stays shLowLevel.tvANSFilter.
        % NOTE: everything goes through ONE 'Name', value cell - MATLAB
        % forbids following name=value syntax with a cell expansion.
        fwd = {'NoiseCov', opts.NoiseCov, ...
               'Constraints', opts.Constraints, ...
               'SignalMode', char(opts.SignalMode), ...
               'NIterSignal', opts.NIterSignal, ...
               'Robust', opts.Robust, ...
               'Blocks', char(opts.Blocks)};
        if ~isempty(opts.Shrinkage)
            fwd = [fwd, {'Shrinkage', opts.Shrinkage}];
        end
        if ~isempty(opts.VCEMinDegree)
            fwd = [fwd, {'VCEMinDegree', opts.VCEMinDegree}];
        end
        if ~isempty(opts.VCEBands)
            fwd = [fwd, {'VCEBands', opts.VCEBands}];
        end
        [Xf, op, info] = shLowLevel.tvANSFilter(X, obj.epochs, idx, fwd{:});
        out = obj;
        out.sigmaCs = nan(size(obj.Cs));
        out.sigmaSs = nan(size(obj.Ss));
        for k = 1:obj.nEpochs
            [cf, sf] = shLowLevel.csFromVec(Xf(:,k), idx);
            [sc, ss] = shLowLevel.csFromVec(info.sigmaXfres(:,k), idx);
            % degrees below MinDegree pass through unfiltered
            out.Cs(:,:,k) = obj.Cs(:,:,k);
            out.Ss(:,:,k) = obj.Ss(:,:,k);
            nz = idx.minDegree + 1;
            out.Cs(nz:end, :, k) = cf(nz:end, :);
            out.Ss(nz:end, :, k) = sf(nz:end, :);
            out.sigmaCs(nz:end, :, k) = sc(nz:end, :);
            out.sigmaSs(nz:end, :, k) = ss(nz:end, :);
        end
        out.history(end+1) = sprintf(...
            "tvANS filtered (SignalMode=%s, constraints=%d)", ...
            opts.SignalMode, size(opts.Constraints, 2));
    end

    function [avg, out] = basinAverage(obj, B, opts)
        %BASINAVERAGE Basin averages, naive or exactly deconvolved.
        %   AVG = ts.basinAverage(B) - naive averages B'*x ./ diag(B'*B),
        %   K x T, with B (P x K) basin kernels in shLowLevel.shIndex ordering
        %   (build from a grid mask via shLowLevel.synthesisMatrix:
        %   b = Y' * (w .* mask)).
        %   [AVG, OUT] = ts.basinAverage(B, Deconvolve (false)=true, Op (struct())=op) removes
        %   filter attenuation and inter-basin leakage exactly using the
        %   stored operator from ts.filter (see shLowLevel.basinDeconvolve).
        %   OUT.sigma (K x T) is the 1-sigma noise uncertainty of AVG
        %   propagated through the deconvolution (v2.1).
        %   Outputs
        %     avg        (K x T) double   basin averages (deconvolved when Deconvolve=true)
        %     out        struct: sigma (K x T), c, attn, condA (deconvolution path)
        %
        %   Options
        %     Ridge (0)  Tikhonov ridge added to the deconvolution normal
        %         matrix; raise it when the kernel matrix is ill-conditioned
        %         (out.condA reports the condition number)
        arguments
            obj
            B double
            opts.Deconvolve (1,1) logical = false
            opts.Op struct = struct()
            opts.Ridge (1,1) double {mustBeNonnegative} = 0   % relative, see shLowLevel.basinDeconvolve
        end
        if opts.Deconvolve
            if isempty(fieldnames(opts.Op))
                error('shSeries:missingOp', ...
                    'Deconvolve=true requires Op=op from ts.filter("tvANS", ...).');
            end
            [avg, out] = shLowLevel.basinDeconvolve(B, opts.Op, Ridge = opts.Ridge);
        else
            idx = shLowLevel.shIndex(obj.nmax);
            if size(B,1) ~= idx.P
                error('shSeries:badBasis', ...
                    'B must be P x K in shLowLevel.shIndex(nmax=%d) ordering (P=%d).', ...
                    obj.nmax, idx.P);
            end
            X = zeros(idx.P, obj.nEpochs);
            for k = 1:obj.nEpochs
                X(:,k) = shLowLevel.vecFromCS(obj.Cs(:,:,k), obj.Ss(:,:,k), idx);
            end
            avg = (B' * X) ./ diag(B' * B);
            out = struct();
        end
    end

    function disp(obj)
        %DISP Compact display: months, span, gaps, processing history.
        %   disp(TS) prints the product type, the number of epochs, the
        %   maximum degree and the covered epoch range, followed by the
        %   processing history. Called automatically when an shSeries
        %   object is shown without a semicolon.
        %
        %   Inputs
        %     obj  (1,1) shSeries  the series to display
        %   Outputs
        %     none - the summary is written to the command window
        %
        %   Example
        %     ts = shSeries.read("ITSG-Grace2018_n60_*.gfc");
        %     disp(ts)
        %
        %   Developed by Matthias Weigelt with the help of Claude (Fable 5).
        fprintf('  shSeries: %s | T=%d epochs | nmax=%d | %.4f..%.4f\n', ...
            obj.productType, obj.nEpochs, obj.nmax, ...
            min(obj.epochs), max(obj.epochs));
        fprintf('    history:\n');
        fprintf('      %s\n', obj.history);
    end
end

methods (Static)
    function obj = read(files, opts)
        %READ Read a series of gfc files (pattern or list), sorted by epoch.
        %   TS = shSeries.read("GSM-2_*.gfc") or shSeries.read(fileList).
        %   Outputs
        %     obj        (1 x 1) shSeries   wildcard-read, epoch-sorted stack
        %
        %   Options
        %     Truncate (NaN)  truncate every field to this nmax while
        %         reading (NaN: keep the file resolution)
        %
        %   Outputs
        %     obj  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            files
            opts.Truncate (1,1) double = NaN
        end
        if (ischar(files) || (isstring(files) && isscalar(files))) ...
                && contains(string(files), ["*", "?"])
            d = dir(char(files));
            if isempty(d)
                error('shSeries:noFiles', 'Pattern matched no files: %s', ...
                    char(files));
            end
            list = string(fullfile({d.folder}, {d.name}))';
            % writeGFC drops "<file>.gfc.provenance.json" next to each
            % export, and the natural pattern "*.gfc*" (which has to be
            % loose enough for ".gfc.gz") matches those too - so a folder
            % written BY the toolbox could not be read back BY the
            % toolbox. Sidecars are metadata, never solutions.
            list = list(~endsWith(list, ".provenance.json"));
            if isempty(list)
                error('shSeries:noFiles', ...
                    'Pattern matched only provenance sidecars: %s', ...
                    char(files));
            end
        else
            list = string(files); list = list(:);
        end
        arr = shCoefficients.empty(0, 1);
        for k = 1:numel(list)
            g = shCoefficients.read(list(k));
            if ~isnan(opts.Truncate), g = g.truncate(opts.Truncate); end
            arr(k, 1) = g;
        end
        obj = shSeries(arr);
    end

    function obj = fromFolder(folder, opts)
        %FROMFOLDER Read every matching gfc file in a folder as a series.
        %   TS = shSeries.fromFolder("itsg_series")
        %   TS = shSeries.fromFolder(dest, Pattern="*n96*.gfc*", Truncate (NaN)=60)
        %   Convenience wrapper around shSeries.read (sorted by epoch);
        %   pairs with shLowLevel.fetchITSG, which downloads ITSG monthly
        %   solutions into tests/test_data/itsg_series on demand.
        %   Outputs
        %     obj        (1 x 1) shSeries   all gfc(.gz) files of a folder, epoch-sorted
        %
        %   Outputs
        %     obj  (1,1) shSeries  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            folder (1,1) string
            opts.Pattern (1,1) string = "*.gfc*"
            opts.Truncate (1,1) double = NaN
        end
        assert(isfolder(folder), 'shSeries:noFiles', ...
            'Not a folder: %s', folder);
        obj = shSeries.read(string(fullfile(char(folder), ...
            char(opts.Pattern))), Truncate = opts.Truncate);
    end
end

methods (Access = private)
    function assertClean(obj, opname)
        if any(isnan(obj.Cs(:))) || any(isnan(obj.Ss(:)))
            error('shSeries:nanInSeries', ...
                '%s requires NaN-free stacks - fill or drop gap epochs first.', ...
                opname);
        end
    end

    function X = flatten(obj)
        Nc = (obj.nmax + 1)^2;
        X = [reshape(obj.Cs, Nc, obj.nEpochs);
             reshape(obj.Ss, Nc, obj.nEpochs)];
    end

    function [Cs, Ss] = unflatten(obj, X)
        Nc = (obj.nmax + 1)^2;
        Cs = reshape(X(1:Nc, :), obj.nmax+1, obj.nmax+1, obj.nEpochs);
        Ss = reshape(X(Nc+1:end, :), obj.nmax+1, obj.nmax+1, obj.nEpochs);
    end
end
end

function [ep, ord] = sortWithNaN(ep)
% sort by epoch, NaNs kept at the end in original order
[~, ord] = sort(ep);
ep = ep(ord);
end
