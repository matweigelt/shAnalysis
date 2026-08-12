classdef shClimatology
%SHCLIMATOLOGY Deterministic climatology of an SH series.
%
%   Bias + trend + annual + semi-annual model per Stokes coefficient,
%   as fitted by shSeries.climatology (least squares or Huber-robust):
%     x(t) = bias + trend*(t-t0) + cosA*cos(2*pi*(t-t0)) + sinA*sin(...)
%                 + cosSA*cos(4*pi*(t-t0)) + sinSA*sin(4*pi*(t-t0))
%                 + sum_k extra cos/sin terms of period p_k  (v2.1)
%   All phases and the bias refer to the reference epoch t0. Extra
%   periods (Periods= option) typically hold the GRACE tidal alias terms
%   S2 (161 d), K2 (3.66 yr), K1 (7.48 yr); access them via periodic(k).
%   Per-coefficient 1-sigma significance from the fit is stored in
%   coefSigma and attached to every component accessor's output, so
%   trend().sigmaC etc. enable significance-masked maps.
%
%   Usage
%     [clim, resid] = ts.climatology(Robust=true);
%     g   = clim.eval(2024.5);          % full model field at an epoch
%     tr  = clim.trend;                 % shCoefficients, per year
%     b0  = clim.bias;                  % field at t0
%     [A, lat, lon] = clim.amplitudeMap("annual", -90:90, 0:359, ...
%                        quantity="ewh", kn=kn);   % pointwise amplitude
%
%   Component accessors return shCoefficients (so the whole toolbox -
%   synthesis, spectra, filtering - applies to them directly).
%
%   Properties (read-only)
%     t0            (1,1) double   reference epoch, decimal years
%     GM, R, nmax   template constants
%     history       string column
%
%   See also shSeries, shCoefficients.
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

properties (SetAccess = private)
    t0 (1,1) double = NaN
    GM (1,1) double = 3.986004415e14
    R (1,1) double = 6378136.3
    biasC double
    biasS double
    trendC double
    trendS double
    cosAnnC double
    cosAnnS double
    sinAnnC double
    sinAnnS double
    cosSemiC double
    cosSemiS double
    sinSemiC double
    sinSemiS double
    periods (1,:) double = double.empty(1,0)   % extra periods [yr] (v2.1)
    extraCosC double = []    % n1 x n1 x K stacks for the extra periods
    extraCosS double = []
    extraSinC double = []
    extraSinS double = []
    coefSigma double = []    % (6+2K) x 2*Nc per-coefficient 1-sigma (v2.1)
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

    function n = get.nmax(obj), n = size(obj.biasC, 1) - 1; end

    function obj = withNote(obj, txt)
        %WITHNOTE Append a provenance note to the history (v2.2).
        %   OBJ = clim.withNote(TXT) - the only supported way to annotate
        %   from outside the class (history is SetAccess=private).
        %   Inputs
        %     txt  free-text provenance note appended to the history
        %
        %   Outputs
        %     obj        (1 x 1) shClimatology   copy with the note appended to its history
        %
        %   Outputs
        %     obj  (1,1) shClimatology  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            txt {mustBeTextScalar}
        end
        obj.history(end+1) = string(txt);
    end

    function out = removeGIA(obj, gia)
        %REMOVEGIA Subtract a GIA rate model from the fitted trend.
        %   OUT = clim.removeGIA(GIA): GIA is an shCoefficients holding
        %   RATE coefficients [1/yr] (user-supplied model file via
        %   shCoefficients.read). Only trendC/trendS change; trend sigmas
        %   are kept - the model is treated as exact (documented
        %   limitation; GIA model spread often dominates regional trend
        %   error budgets - compare several models).
        %   Inputs
        %     gia  GIA rate field (shCoefficients, struct or file) removed
        %          from the trend component
        %
        %   Outputs
        %     out  (1,1) shClimatology  copy with the GIA trend removed from the
        %          trend component; history appended
        %
        %   Outputs
        %     out  (1,1) shClimatology  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            obj
            gia (1,1) shCoefficients
        end
        n1 = obj.nmax + 1;
        nc = min(n1, size(gia.C, 1));
        out = obj;
        out.trendC(1:nc, 1:nc) = obj.trendC(1:nc, 1:nc) - gia.C(1:nc, 1:nc);
        out.trendS(1:nc, 1:nc) = obj.trendS(1:nc, 1:nc) - gia.S(1:nc, 1:nc);
        out.history(end+1) = sprintf( ...
            "GIA removed from trend (model nmax=%d, treated exact)", nc - 1);
    end

    function g = eval(obj, epoch)
        %EVAL Evaluate the full climatology at an epoch -> shCoefficients.
        %   Outputs
        %     g          (1 x 1) shCoefficients   model field evaluated at the epoch
        %
        %   Outputs
        %     g  (1,1) shCoefficients  climatology evaluated at the epoch
        arguments
            obj
            epoch (1,1) double
        end
        dt = epoch - obj.t0;
        ca = cos(2*pi*dt); sa = sin(2*pi*dt);
        cs = cos(4*pi*dt); ss = sin(4*pi*dt);
        C = obj.biasC + obj.trendC*dt + obj.cosAnnC*ca + obj.sinAnnC*sa ...
            + obj.cosSemiC*cs + obj.sinSemiC*ss;
        S = obj.biasS + obj.trendS*dt + obj.cosAnnS*ca + obj.sinAnnS*sa ...
            + obj.cosSemiS*cs + obj.sinSemiS*ss;
        for k = 1:numel(obj.periods)
            cp = cos(2*pi*dt/obj.periods(k));
            sp = sin(2*pi*dt/obj.periods(k));
            C = C + obj.extraCosC(:,:,k)*cp + obj.extraSinC(:,:,k)*sp;
            S = S + obj.extraCosS(:,:,k)*cp + obj.extraSinS(:,:,k)*sp;
        end
        g = shCoefficients(C, S, GM = obj.GM, R = obj.R, Epoch = epoch, ...
            ProductType = "climatology", Name = "climatology", ...
            History = [obj.history; sprintf("evaluated at %.4f", epoch)]);
    end

    function g = bias(obj)
        %BIAS Field at the reference epoch t0.
        %   Outputs
        %     g          (1 x 1) shCoefficients   bias component with sigmas
        %
        %   Outputs
        %     out  (1,1) shCoefficients  bias (mean) component; sigmas from the fit
        g = obj.component(obj.biasC, obj.biasS, "climatology bias (at t0)");
    end

    function g = trend(obj)
        %TREND Linear trend field, units of the coefficients per year.
        %   Outputs
        %     g          (1 x 1) shCoefficients   trend component [units/yr] with sigmas
        %
        %   Outputs
        %     out  (1,1) shCoefficients  trend component [1/yr]; sigmas from the fit
        g = obj.component(obj.trendC, obj.trendS, "climatology trend [1/yr]");
    end

    function g = cosAnnual(obj)
        %COSANNUAL Cosine annual component (sigmas attached).
        %   Outputs
        %     g          (1 x 1) shCoefficients   cosine annual component
        %
        %   Outputs
        %     out  (1,1) shCoefficients  component field; sigmas propagated from the fit
        g = obj.component(obj.cosAnnC, obj.cosAnnS, "annual cos component");
    end
    function g = sinAnnual(obj)
        %SINANNUAL Sine annual component (sigmas attached).
        %   Outputs
        %     g          (1 x 1) shCoefficients   sine annual component
        %
        %   Outputs
        %     out  (1,1) shCoefficients  component field; sigmas propagated from the fit
        g = obj.component(obj.sinAnnC, obj.sinAnnS, "annual sin component");
    end
    function g = cosSemiannual(obj)
        %COSSEMIANNUAL Cosine semiannual component (sigmas attached).
        %   Outputs
        %     g          (1 x 1) shCoefficients   cosine semiannual component
        %
        %   Outputs
        %     out  (1,1) shCoefficients  component field; sigmas propagated from the fit
        g = obj.component(obj.cosSemiC, obj.cosSemiS, "semi-annual cos component");
    end
    function g = sinSemiannual(obj)
        %SINSEMIANNUAL Sine semiannual component (sigmas attached).
        %   Outputs
        %     g          (1 x 1) shCoefficients   sine semiannual component
        %
        %   Outputs
        %     out  (1,1) shCoefficients  component field; sigmas propagated from the fit
        g = obj.component(obj.sinSemiC, obj.sinSemiS, "semi-annual sin component");
    end

    function [gc, gs] = periodic(obj, k)
        %PERIODIC Cos/sin components of the k-th extra period (v2.1).
        %   [GC, GS] = clim.periodic(K) returns shCoefficients for the
        %   cos and sin terms of period obj.periods(K) [yr] - e.g. the
        %   S2/K2/K1 tidal alias terms when fitted via
        %   ts.climatology(Periods=[161/365.25, 3.66, 7.48]).
        %   Inputs
        %     k  harmonic index: 1 = annual, 2 = semiannual, ...
        %
        %   Outputs
        %     gc         (1 x 1) shCoefficients   cosine component of the k-th extra period
        %     gs         (1 x 1) shCoefficients   sine component
        %
        %   Outputs
        %     out  (1,1) shCoefficients  cos/sin component of the k-th extra period
        arguments
            obj
            k (1,1) double {mustBeInteger, mustBePositive}
        end
        if k > numel(obj.periods)
            error('shClimatology:badIndex', ...
                'Only %d extra periods are stored.', numel(obj.periods));
        end
        lbl = sprintf("period %.6g yr", obj.periods(k));
        gc = obj.component(obj.extraCosC(:,:,k), obj.extraCosS(:,:,k), ...
            lbl + " cos", 5 + 2*k);
        gs = obj.component(obj.extraSinC(:,:,k), obj.extraSinS(:,:,k), ...
            lbl + " sin", 6 + 2*k);
    end

    function [A, lat, lon, phase] = amplitudeMap(obj, which, latVec, lonVec, opts)
        %AMPLITUDEMAP Pointwise seasonal amplitude (and phase) map.
        %   [A, LAT, LON, PHASE] = clim.amplitudeMap(WHICH, ... with WHICH
        %   "annual" | "semiannual" | k (index into clim.periods, v2.1),
        %   LATVEC, LONVEC, quantity ("geoid")=..., kn ([])=...) synthesizes the cos and
        %   sin component fields and combines them pointwise:
        %       A     = sqrt(fc.^2 + fs.^2)
        %       PHASE = atan2(fs, fc)   [rad, relative to t0]
        %   (amplitude of a harmonic is NOT a linear functional of the
        %   coefficients, so it must be formed in the space domain).
        %   Outputs
        %     A          (nlat x nlon) double   amplitude of the harmonic in the requested quantity
        %     lat        (1 x nlat) double   latitudes
        %     lon        (1 x nlon) double   longitudes
        %     phase      (nlat x nlon) double   phase [rad] of the harmonic
        %
        %   Outputs
        %     A  (nlat x nlon) double  annual amplitude of the synthesized quantity
        %     h  (1,1) graphics handle  figure (Plot = true only)
        arguments
            obj
            which (1,1)   % "annual" | "semiannual" | integer k (extra period)
            latVec (1,:) double
            lonVec (1,:) double
            opts.quantity (1,1) string = "geoid"
            opts.kn double = []
        end
        if isnumeric(which)
            [gc, gs] = obj.periodic(which);
        elseif string(which) == "annual"
            gc = obj.cosAnnual; gs = obj.sinAnnual;
        elseif string(which) == "semiannual"
            gc = obj.cosSemiannual; gs = obj.sinSemiannual;
        else
            error('shClimatology:badInput', ...
                'which must be "annual", "semiannual" or an extra-period index.');
        end
        args = {latVec, lonVec, 'quantity', opts.quantity};
        if ~isempty(opts.kn), args = [args, {'kn', opts.kn}]; end
        [fc, lat, lon] = gc.synthesis(args{:});
        fs = gs.synthesis(args{:});
        A = sqrt(fc.^2 + fs.^2);
        phase = atan2(fs, fc);
    end

    function disp(obj)
        %DISP Compact display: fitted components, reference epoch, note.
        %   disp(CLIM) prints the maximum degree and the reference epoch
        %   t0, the list of fitted components (bias, trend, annual and
        %   semi-annual cosine/sine, plus any extra Periods) and the
        %   processing history including the GIA note. Called
        %   automatically when an shClimatology object is shown without
        %   a semicolon.
        %
        %   Inputs
        %     obj  (1,1) shClimatology  the fitted model to display
        %   Outputs
        %     none - the summary is written to the command window
        %
        %   Example
        %     clim = ts.climatology();
        %     disp(clim)
        %
        %   Developed by Matthias Weigelt with the help of Claude (Fable 5).
        fprintf('  shClimatology: nmax=%d | t0=%.4f\n', obj.nmax, obj.t0);
        fprintf('    components: bias, trend, annual (cos/sin), semi-annual (cos/sin)\n');
        fprintf('    history:\n');
        fprintf('      %s\n', obj.history);
    end
end

methods (Static)
    function obj = fromCoef(coef, t0, series, opts)
        %FROMCOEF Build from shLowLevel.fitDeterministicModel output (internal).
        %   OBJ = shClimatology.fromCoef(COEF, T0, SERIES, Periods (double.empty(1,0))=[],
        %   CoefSigma=[]) with COEF (6+2K) x (2*Nc) over the flattened
        %   [C(:); S(:)] stack of SERIES; rows 7:end are the cos/sin
        %   pairs of the K extra Periods [yr].
        %   Inputs
    %     coef   (P x K) fitted component coefficients in shIndex ordering
    %     t0     (1 x 1) reference epoch of the fit [decimal years]
    %     series (1 x 1) shSeries the fit came from (history carrier)
    %   Options
    %     Periods ([])      (1 x K2) extra harmonic periods [yr]
    %     CoefSigma ([])    (P x K) 1-sigma uncertainties of coef
    %
    %   Outputs
        %     obj        (1 x 1) shClimatology   rebuilt from a coefficient table with the series' layout
        %
        %   Outputs
        %     obj  (1,1) shClimatology  modified copy; the operation is appended to the
        %              history (immutable value-class pattern)
        arguments
            coef double
            t0 (1,1) double
            series (1,1) shSeries
            opts.Periods (1,:) double = double.empty(1,0)
            opts.CoefSigma double = []
        end
        K = numel(opts.Periods);
        assert(size(coef,1) == 6 + 2*K, 'shClimatology:badInput', ...
            'coef must have %d rows for %d extra periods.', 6+2*K, K);
        n1 = series.nmax + 1;
        Nc = n1^2;
        obj = shClimatology();
        obj.t0 = t0; obj.GM = series.GM; obj.R = series.R;
        names = {'bias','trend','cosAnn','sinAnn','cosSemi','sinSemi'};
        for r = 1:6
            obj.([names{r} 'C']) = reshape(coef(r, 1:Nc),      n1, n1);
            obj.([names{r} 'S']) = reshape(coef(r, Nc+1:end),  n1, n1);
        end
        obj.periods = opts.Periods;
        if K > 0
            obj.extraCosC = zeros(n1, n1, K); obj.extraCosS = zeros(n1, n1, K);
            obj.extraSinC = zeros(n1, n1, K); obj.extraSinS = zeros(n1, n1, K);
            for k = 1:K
                rc = 5 + 2*k; rs = 6 + 2*k;
                obj.extraCosC(:,:,k) = reshape(coef(rc, 1:Nc),     n1, n1);
                obj.extraCosS(:,:,k) = reshape(coef(rc, Nc+1:end), n1, n1);
                obj.extraSinC(:,:,k) = reshape(coef(rs, 1:Nc),     n1, n1);
                obj.extraSinS(:,:,k) = reshape(coef(rs, Nc+1:end), n1, n1);
            end
        end
        obj.coefSigma = opts.CoefSigma;
        obj.history = [series.history; ...
            sprintf("climatology fitted (t0=%.4f, T=%d epochs, %d extra periods)", ...
                t0, series.nEpochs, K)];
    end
end

methods (Access = private)
    function g = component(obj, C, S, label, row)
        if nargin < 5, row = obj.componentRow(label); end
        args = {'GM', obj.GM, 'R', obj.R, 'Epoch', obj.t0, ...
            'ProductType', "climatology", 'Name', label, ...
            'History', [obj.history; label]};
        if ~isempty(obj.coefSigma) && row > 0
            n1 = obj.nmax + 1; Nc = n1^2;
            sC = reshape(obj.coefSigma(row, 1:Nc),     n1, n1);
            sS = reshape(obj.coefSigma(row, Nc+1:end), n1, n1);
            args = [args, {'SigmaC', sC, 'SigmaS', sS}];
        end
        g = shCoefficients(C, S, args{:});
    end

    function r = componentRow(obj, label)
        % map accessor labels to coef rows (0 = unknown)
        lbl = char(label);
        if     contains(lbl, 'bias'),            r = 1;
        elseif contains(lbl, 'trend'),           r = 2;
        elseif contains(lbl, 'semi-annual cos'), r = 5;
        elseif contains(lbl, 'semi-annual sin'), r = 6;
        elseif contains(lbl, 'annual cos'),      r = 3;
        elseif contains(lbl, 'annual sin'),      r = 4;
        else, r = 0;
        end
    end
end
end
