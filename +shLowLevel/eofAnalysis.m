function [modes, pcs, varExplained, info] = eofAnalysis(ts, opts)
%EOFANALYSIS Empirical orthogonal functions of a coefficient series.
%
%   [MODES, PCS, VAREXPLAINED, INFO] = shLowLevel.eofAnalysis(TS) decomposes
%   the (typically residual) series into spatial modes and principal
%   components via SVD of the flattened, time-centered coefficient
%   stack: X = U S V'; MODES{k} carries U(:,k)*S(k,k)/sqrt(T-1) as an
%   shCoefficients (synthesize for the spatial pattern), PCS(:,k) =
%   sqrt(T-1)*V(:,k) is the unit-variance temporal amplitude, so
%   mode-k field * pc reconstructs the mode's contribution exactly.
%
%   Standard use: identify leading interannual patterns (ENSO-like
%   hydrology, ice-sheet acceleration) after removing the climatology:
%       [~, resid] = ts.climatology();
%       [modes, pcs, ve] = shLowLevel.eofAnalysis(resid);
%
%   Inputs
%     ts   (1,1) shSeries (NaN-free; use dropNaN across the gap)
%   Options
%     NModes (6)         modes returned
%     Center (true)      remove the temporal mean first
%   Outputs
%     modes         1 x K cell of shCoefficients (spatial patterns)
%     pcs           (T,K) double, unit variance columns
%     varExplained  (K,1) double, fraction in [0,1]
%     info          struct: singular values, total variance, T
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    ts (1,1) shSeries
    opts.NModes (1,1) double {mustBeInteger, mustBePositive} = 6
    opts.Center (1,1) logical = true
end
assert(~any(isnan(ts.Cs(:))) && ~any(isnan(ts.Ss(:))), ...
    'shLowLevel:eofAnalysis:nanInSeries', ...
    'EOF analysis needs a NaN-free series - bridge the gap with dropNaN.');
n1 = ts.nmax + 1; Nc = n1^2; T = ts.nEpochs;
X = [reshape(ts.Cs, Nc, T); reshape(ts.Ss, Nc, T)];  % (2Nc) x T
assert(T >= 3, 'shLowLevel:eofAnalysis:tooFewEpochs', ...
    'EOF analysis needs >= 3 epochs (got %d).', T);
if opts.Center
    X = X - mean(X, 2);
end
[U, Sv, V] = svd(X, 'econ');
sv = diag(Sv);
K = min(opts.NModes, numel(sv));
totVar = sum(sv.^2);
varExplained = sv(1:K).^2 / max(totVar, realmin);
modes = cell(1, K);
for k = 1:K
    u = U(:, k) * sv(k) / sqrt(max(T - 1, 1));
    Ck = reshape(u(1:Nc), n1, n1);
    Sk = reshape(u(Nc+1:end), n1, n1);
    modes{k} = shCoefficients(Ck, Sk, GM = ts.GM, R = ts.R, ...
        Name = sprintf("EOF mode %d (%.1f%%)", k, 100*varExplained(k)));
end
pcs = sqrt(max(T - 1, 1)) * V(:, 1:K);
info = struct('singularValues', sv, 'totalVariance', totVar, 'T', T, ...
    'centered', opts.Center);
end
