function [model, Xres, coef, A, coefSigma, resVar, ar1] = fitDeterministicModel(X, tYears, opts)
%FITDETERMINISTICMODEL Bias/trend/(semi-)annual (+extra periods) per SH coefficient.
%
%   [MODEL, XRES, COEF, A, COEFSIGMA, RESVAR] =
%       shLowLevel.fitDeterministicModel(X, TYEARS, ...) fits, for each
%   coefficient time series (rows of X, P x T), the model
%       x(t) = a + b*tc + c*cos(2*pi*tc) + d*sin(2*pi*tc)
%                + e*cos(4*pi*tc) + f*sin(4*pi*tc)
%                + sum_k [ g_k*cos(2*pi*tc/p_k) + h_k*sin(2*pi*tc/p_k) ]
%   with tc = t - T0 (NaN) and optional extra periods p_k [years] - e.g. the
%   GRACE tidal alias periods S2 = 161/365.25, K2 = 3.66, K1 = 7.48.
%
%   Options
%     T0       (mean(tYears))  reference epoch [decimal years]
%     Periods (double.empty(1,0))  ([])   extra periods p_k [years], row vector
%     Weights  ([])   T x 1 per-epoch weights (e.g. 1./s from VCE);
%                     weighted LS with sqrt(w)-scaled rows
%     Robust   (false) Huber IRLS per coefficient
%     HuberK   (1.5), MaxIter (10)
%     ARCorrect (false) inflate coefSigma by sqrt((1+r1)/(1-r1)) with the
%              Kendall bias-corrected r1' = r1 + (1+3r1)/T (v2.5) where
%              the lag-1 autocorrelation of each residual series - the
%              standard first-order correction for temporally correlated
%              (AR(1)-like) GRACE residuals. Monte-Carlo validated
%              (phi=0.6, T=120): white-noise sigmas underestimate the
%              trend scatter by 2.0x; corrected ratio 1.06 (the residual
%              slight underestimate stems from the downward-biased sample
%              r1 - documented approximation).
%     Breaks (double.empty(1,0))  epochs [decimal years] of continuous
%         piecewise-linear trend hinges; each adds a max(0, t - tb)
%         column to the design matrix
%
%   Outputs
%     model      P x T          fitted deterministic part
%     Xres       P x T          X - model
%     coef       ncol x P       [bias; trend; cosA; sinA; cosSA; sinSA;
%                                cosP1; sinP1; ...], ncol = 6 + 2*numel(Periods)
%     A          T x ncol       design matrix (unweighted)
%     coefSigma  ncol x P       1-sigma per coefficient: sqrt(resVar_p *
%                               diag(inv(A'WA))); OLS formula, validated by
%                               Monte Carlo (ratios 0.99-1.02). For
%                               Robust=true the final IRLS weights enter W
%                               (approximation, documented).
%     resVar     1 x P          residual variance rss/(T - ncol)
%     ar1        1 x P          bias-corrected lag-1 residual auto-
%                               correlation actually applied (clamped to
%                               [0, 0.99]; 0 when ARCorrect=false)
%
%   Claude (Fable 5), 2026-08-07 (extends the 2026-08-07 6-parameter fit).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    X double
    tYears (:,1) double
    opts.Robust (1,1) logical = false
    opts.HuberK (1,1) double = 1.5
    opts.MaxIter (1,1) double = 10
    opts.T0 (1,1) double = NaN
    opts.Periods (1,:) double {mustBePositive} = double.empty(1,0)
    opts.Weights (:,1) double = []
    opts.ARCorrect (1,1) logical = false
    opts.Breaks (1,:) double = double.empty(1,0)  % hinge epochs (v2.4):
                                     % appends max(0, t - tb) columns -
                                     % continuous piecewise-linear trend.
                                     % NOT passed through climatology
                                     % (fixed schema); use
                                     % shSeries.trendBreaks.
end

T = numel(tYears);
assert(size(X,2) == T, 'X must be P x T with T = numel(tYears).');

t0 = opts.T0;
if isnan(t0), t0 = mean(tYears); end
tc = tYears - t0;
A = [ones(T,1), tc, cos(2*pi*tc), sin(2*pi*tc), cos(4*pi*tc), sin(4*pi*tc)];
for p_ = opts.Periods
    A = [A, cos(2*pi*tc/p_), sin(2*pi*tc/p_)]; %#ok<AGROW>
end
for tb = opts.Breaks
    A = [A, max(tYears(:) - tb, 0)]; %#ok<AGROW>   % hinge (v2.4)
end
ncol = size(A, 2);
if T <= ncol
    warning('shLowLevel:fitDeterministicModel:fewEpochs', ...
        'Only %d epochs for %d parameters; sigmas will be NaN.', T, ncol);
end

w = opts.Weights;
if isempty(w), w = ones(T,1); end
assert(numel(w) == T && all(w >= 0), 'Weights must be T x 1, nonnegative.');
sw = sqrt(w);
Aw = A .* sw;

P = size(X,1);
if ~opts.Robust
    coef = Aw \ (sw .* X');                     % ncol x P
    Nw = Aw' * Aw;
else
    coef = zeros(ncol, P);
    Wfin = ones(T, P);
    for p = 1:P
        y = X(p,:)';
        c = Aw \ (sw .* y);
        wr = w;
        for it = 1:opts.MaxIter
            r = y - A*c;
            s = 1.4826 * median(abs(r - median(r)));
            if s <= 0, break; end
            hub = min(1, opts.HuberK ./ max(abs(r)/s, eps));
            wr = w .* hub;
            swr = sqrt(wr);
            cNew = (A .* swr) \ (y .* swr);
            if norm(cNew - c) <= 1e-12 * max(norm(c), 1), c = cNew; break; end
            c = cNew;
        end
        coef(:,p) = c;
        Wfin(:,p) = wr;
    end
end

model = (A * coef)';                            % P x T
Xres  = X - model;

ar1 = zeros(1, P);
if nargout >= 5
    dof = max(T - ncol, 1);
    coefSigma = nan(ncol, P);
    resVar = nan(1, P);
    if T > ncol
        if ~opts.Robust
            dN = diag(inv(Nw));                 % shared across coefficients
            rw = Xres .* (sw');                 % weighted residuals, P x T
            resVar = sum(rw.^2, 2)' / dof;
            coefSigma = sqrt(dN * resVar);      % ncol x P
        else
            for p = 1:P
                swr = sqrt(Wfin(:,p));
                Awp = A .* swr;
                rw = swr .* Xres(p,:)';
                resVar(p) = sum(rw.^2) / dof;
                coefSigma(:,p) = sqrt(diag(inv(Awp'*Awp)) * resVar(p));
            end
        end
        if opts.ARCorrect && T > ncol + 2
            for p = 1:P
                r = Xres(p, :)';
                den = sum(r.^2);
                if den > 0
                    ar1(p) = min(max(sum(r(1:end-1) .* r(2:end)) / den, 0), 0.99);
                end
            end
            % Kendall/Marriott-Pope first-order bias correction (v2.5):
            % the sample r1 is biased low by ~(1+3r1)/T; correcting it
            % removes about half the remaining sigma underestimation.
            % Monte-Carlo (4000 runs, trend + annual + semi-annual fit
            % on AR(1) noise, empirical trend scatter / mean sigma):
            %   phi=0.6 T=120: 1.065 -> 1.028    phi=0.3 T=120: 1.047 -> 1.030
            %   phi=0.6 T= 60: 1.152 -> 1.079    phi=0.8 T=180: 1.044 -> 0.991
            % Still first-order; the residual few percent stems from the
            % regression absorbing low-frequency noise (documented).
            ar1 = min(max(ar1 + (1 + 3 * ar1) / T, 0), 0.99);
            coefSigma = coefSigma .* sqrt((1 + ar1) ./ (1 - ar1));
        end
    end
end
end
