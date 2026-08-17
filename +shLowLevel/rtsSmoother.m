function smo = rtsSmoother(filt)
%RTSSMOOTHER Rauch-Tung-Striebel backward pass for kalmanFilter output.
%
%   SMO = shLowLevel.rtsSmoother(FILT) runs the fixed-interval RTS
%   smoother (Kurtenbach 2012, Sec. 3.3.3) over the forward-filter
%   result of shLowLevel.kalmanFilter:
%       G_t  = Pf_t B' Pp_{t+1}^-1
%       xs_t = xf_t + G_t (xs_{t+1} - xp_{t+1})
%       Ps_t = Pf_t + G_t (Ps_{t+1} - Pp_{t+1}) G_t'
%   so the estimate at epoch t uses ALL observations, past and future.
%   With the stationary initialization of kalmanFilter this reproduces
%   the joint least-squares adjustment over all epochs exactly (Kvas
%   2019, Sec. 2.3; unit-tested to machine precision), i.e. the
%   forward-backward pass is an exact solver for the block-tridiagonal
%   joint normal-equation system - never build that system explicitly.
%
%   Inputs
%     filt  (1 x 1) struct  from shLowLevel.kalmanFilter with
%           StoreCov = "full" (the smoother needs the full predicted
%           and filtered covariances; "diag" runs error loudly)
%
%   Outputs
%     smo  (1 x 1) struct  smoothed results, physical block only:
%          xs (P x T double) smoothed states, sig (P x T double)
%          formal 1-sigma per coefficient (sqrt of the smoothed
%          covariance diagonal), PsLast (P x P double) full smoothed
%          covariance of the LAST epoch (equals the filtered one
%          there), P (1 x 1), order (1 x 1)
%
%   Example
%     filt = shLowLevel.kalmanFilter(model, obs);   % StoreCov="full"
%     smo  = shLowLevel.rtsSmoother(filt);
%     plot(ep, smo.xs(1, :))                        % first coefficient
%
%   Reference: Rauch, Tung & Striebel (1965); Kurtenbach, DGK C-683
%   (2012), Sec. 3.3.3; Kvas, TU Graz PhD thesis (2019), Sec. 2.3.
%   Numerics pre-validated in Python (tools/dev/validate_kalman.py).
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-17, 20:10 UTC.

arguments
    filt (1,1) struct
end

if filt.storeCov ~= "full"
    error('shLowLevel:rtsSmoother:needFullCov', ...
        ['rtsSmoother needs kalmanFilter(..., StoreCov="full"); the ' ...
         'backward gain uses the full predicted covariances.']);
end

B = filt.B;
P = filt.P;
Pc = size(filt.xf, 1);
T = size(filt.xf, 2);

xs = filt.xf;                          % companion space, overwritten backwards
Ps = filt.Pf(:, :, T);
sig = zeros(P, T);
sig(:, T) = sqrt(max(diag(Ps(1:P, 1:P)), 0));
xsNext = xs(:, T);
PsNext = Ps;
for t = T-1:-1:1
    PfT = filt.Pf(:, :, t);
    PpN = filt.Pp(:, :, t+1);
    G = (PfT * B.') / PpN;             % backward gain, no explicit inverse
    xsT = filt.xf(:, t) + G * (xsNext - filt.xp(:, t+1));
    PsT = PfT + G * (PsNext - PpN) * G.';
    PsT = (PsT + PsT.') / 2;
    xs(:, t) = xsT;
    sig(:, t) = sqrt(max(diag(PsT(1:P, 1:P)), 0));
    xsNext = xsT; PsNext = PsT;
end

smo = struct('xs', xs(1:P, :), 'sig', sig, ...
    'PsLast', filt.Pf(1:P, 1:P, T), 'P', P, 'order', filt.order);
end
