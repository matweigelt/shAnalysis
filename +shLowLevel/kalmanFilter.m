function filt = kalmanFilter(model, obs, opts)
%KALMANFILTER Forward Kalman filter for SH series (Kurtenbach/Kvas).
%
%   FILT = shLowLevel.kalmanFilter(MODEL, OBS) runs the forward Kalman
%   filter of Kurtenbach (2012, Sec. 3.3) with the VAR(p) process model
%   from shLowLevel.estimateVAR and one observation record per epoch.
%   Two observation modes, decided per epoch by the fields present:
%     solution  l = x + v, v ~ N(0, R): a (noisy) coefficient-vector
%               solution with its error covariance - e.g. a GFC series
%               with formal sigmas, or a SINEX COVA solution
%     neq       N x = b: an unreduced normal equation - the SINEX NEQ
%               blocks shLowLevel.readSINEX returns (kind = 'NEQ');
%               applied as information-form update, no R needed
%   An empty record is a data gap: prediction only, the state relaxes
%   toward the process prior (zero mean, covariance -> Sigma0).
%
%   Initialization follows Kvas (2019), eqs. (2.90)-(2.92): x0- = 0,
%   P0- = Sigma0 (the stationary process covariance). With this choice
%   the filter followed by shLowLevel.rtsSmoother is IDENTICAL to the
%   joint least-squares adjustment over all epochs (his Sec. 2.3) -
%   unit-tested to machine precision. A single update with Phi = 0 is
%   the per-epoch Wiener filter, i.e. the temporal generalization of
%   tvANSFilter (also unit-tested).
%
%   Inputs
%     model  (1 x 1) struct  from shLowLevel.estimateVAR: Phi, Q,
%            Sigma0, order, P
%     obs    (1 x T) struct  per-epoch records with fields
%            l (P x 1 double) and R (P x P double, or P x 1 diagonal)
%            for solution mode, or N (P x P double) and b (P x 1
%            double) for neq mode; leave all fields [] for a gap
%
%   Options
%     StoreCov ("full")  "full" keeps the predicted and filtered
%            covariances (needed by rtsSmoother; memory: 2 T (pP)^2
%            doubles - n40, p=1, T=365: ~16 GB; keep T or nmax small,
%            or filter in windows) | "diag" keeps only diagonals
%            (filter-only runs, e.g. NRT one-pass a la Kvas Ch. 3)
%     Contribution (true)  record diag(K) (solution) / diag(P+ N)
%            (neq): the share of the estimate coming from the DATA
%            rather than the process model - Kurtenbach Sec. 3.3.2
%
%   Outputs
%     filt  (1 x 1) struct  companion-space filter results:
%           xf, xp (pP x T double) filtered / predicted states,
%           Pf, Pp (pP x pP x T double, StoreCov="full") or
%           dPf, dPp (pP x T double, StoreCov="diag") covariances,
%           contrib (P x T double) data share per coefficient (NaN in
%           gaps), gap (1 x T logical), P (1 x 1) block size, order
%           (1 x 1), storeCov (1 x 1) string. The physical state is
%           the first block: xf(1:P, :).
%
%   Example
%     model = shLowLevel.estimateVAR(X, Order=1);
%     obs = repmat(struct('l',[],'R',[],'N',[],'b',[]), 1, numel(ep));
%     for t = 1:numel(ep)          % fill l/R (or N/b), leave gaps empty
%         obs(t).l = L(:, t);  obs(t).R = sig(:, t).^2;
%     end
%     filt = shLowLevel.kalmanFilter(model, obs);
%
%   Reference: Kurtenbach, DGK C-683 (2012), Sec. 3.3; Kvas, TU Graz
%   PhD thesis (2019), Secs. 2.2-2.3.
%   Numerics pre-validated in Python (tools/dev/validate_kalman.py).
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-17, 20:10 UTC.

arguments
    model (1,1) struct
    obs (1,:) struct
    opts.StoreCov (1,1) string {mustBeMember(opts.StoreCov, ["full","diag"])} = "full"
    opts.Contribution (1,1) logical = true
end

P = model.P;
p = model.order;
Pc = p * P;
T = numel(obs);

% ---- companion form (Kvas Sec. 2.3.1)
B = zeros(Pc);
B(1:P, :) = [model.Phi{:}];
if p > 1
    B(P+1:end, 1:(p-1)*P) = eye((p-1)*P);
end
Q = zeros(Pc);
Q(1:P, 1:P) = model.Q;
S0 = kron(eye(p), model.Sigma0);        % block-diagonal stationary init
% (exact for p = 1; for p > 1 the block-diagonal Sigma0 approximates the
% companion stationary covariance - conservative, converges in few steps)

full = opts.StoreCov == "full";
xf = zeros(Pc, T); xp = zeros(Pc, T);
if full
    Pf = zeros(Pc, Pc, T); Pp = zeros(Pc, Pc, T);
else
    dPf = zeros(Pc, T); dPp = zeros(Pc, T);
end
contrib = NaN(P, T);
gap = false(1, T);

xPrev = zeros(Pc, 1);
PPrev = S0;
for t = 1:T
    % ---- predict
    xm = B * xPrev;
    Pm = B * PPrev * B.' + Q;
    Pm = (Pm + Pm.') / 2;
    xp(:, t) = xm;
    if full, Pp(:, :, t) = Pm; else, dPp(:, t) = diag(Pm); end

    % ---- update
    o = obs(t);
    hasSol = isfield(o, 'l') && ~isempty(o.l);
    hasNeq = isfield(o, 'N') && ~isempty(o.N);
    if hasSol && hasNeq
        error('shLowLevel:kalmanFilter:ambiguousObs', ...
            'Epoch %d carries both l/R and N/b - provide one mode per epoch.', t);
    end
    if hasSol
        l = o.l(:);
        R = o.R;
        if isvector(R), R = diag(R); end
        if numel(l) ~= P || ~isequal(size(R), [P P])
            error('shLowLevel:kalmanFilter:sizeMismatch', ...
                'Epoch %d: l must be P x 1 and R P x P (or P x 1), P = %d.', t, P);
        end
        % observation touches the first block only: H = [I 0 ... 0]
        PmH = Pm(:, 1:P);                              % Pm * H'
        K = PmH / ((Pm(1:P, 1:P) + R).');              % gain, no inverse
        xPl = xm + K * (l - xm(1:P));
        PPl = Pm - K * PmH.';
        if opts.Contribution
            contrib(:, t) = diag(K(1:P, :));
        end
    elseif hasNeq
        if ~isequal(size(o.N), [P P]) || numel(o.b) ~= P
            error('shLowLevel:kalmanFilter:sizeMismatch', ...
                'Epoch %d: N must be P x P and b P x 1, P = %d.', t, P);
        end
        Nc = zeros(Pc); Nc(1:P, 1:P) = o.N;            % embed in companion space
        bc = zeros(Pc, 1); bc(1:P) = o.b(:);
        Jm = inv(Pm);                                  % information form
        PPl = inv(Jm + Nc);
        xPl = PPl * (Jm * xm + bc); %#ok<MINV>
        if opts.Contribution
            contrib(:, t) = diag(PPl(1:P, 1:P) * o.N);
        end
    else
        gap(t) = true;                                 % prediction only
        xPl = xm; PPl = Pm;
    end
    PPl = (PPl + PPl.') / 2;
    xf(:, t) = xPl;
    if full, Pf(:, :, t) = PPl; else, dPf(:, t) = diag(PPl); end
    xPrev = xPl; PPrev = PPl;
end

filt = struct('xf', xf, 'xp', xp, 'contrib', contrib, 'gap', gap, ...
    'P', P, 'order', p, 'storeCov', opts.StoreCov, 'B', B);
if full
    filt.Pf = Pf; filt.Pp = Pp;
else
    filt.dPf = dPf; filt.dPp = dPp;
end
end
