function out = threeCorneredHat(X)
%THREECORNEREDHAT Per-solution noise levels from N >= 3 series.
%
%   OUT = shx.threeCorneredHat(X) with X (T x N, N >= 3) estimates the
%   individual noise standard deviation of each series by the classic
%   (generalized) three-cornered hat: pairwise difference variances
%   V(i,j) = var(x_i - x_j) cancel the common signal, and under the
%   assumption of mutually independent noises V(i,j) = s_i^2 + s_j^2;
%   the N unknowns are solved from the N*(N-1)/2 pairs by least squares
%   (exact Grubbs solution for N = 3). This answers the question
%   pairwise metrics cannot: WHICH solution is noisy. Negative variance
%   estimates (possible for short series or correlated noises) are
%   clipped to zero and flagged. Python-validated: recovers set noise
%   levels [1, 2, 3, 1.5] within 1% at T = 20000.
%
%   Outputs
%     out        (1,1) struct  fields:
%                  .sigma    (1,N) double  noise standard deviation per series
%                  .variance (1,N) double  raw LS variance estimates (may be < 0)
%                  .pairVar  (P,3) double  [i j var(x_i - x_j)] per pair
%                  .clipped  (1,N) logical true where variance was clipped to 0
%
%   Example
%     out = shx.threeCorneredHat(Y3);   % T x 3 basin series (CSR/GFZ/JPL)
%     fprintf("noise levels: %s\n", sprintf("%.3g ", out.sigma))
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.6.0).
arguments
    X (:,:) double
end
[T, N] = size(X);
if N < 3
    error('shx:threeCorneredHat:needThree', ...
        'Need at least 3 series (got %d).', N);
end
if T < 12
    error('shx:threeCorneredHat:tooShort', ...
        'Need at least 12 epochs (got %d).', T);
end
nP = N * (N - 1) / 2;
A = zeros(nP, N); V = zeros(nP, 1); pairVar = zeros(nP, 3);
k = 0;
for i = 1:N-1
    for j = i+1:N
        k = k + 1;
        d = X(:, i) - X(:, j);
        d = d(isfinite(d));
        V(k) = var(d);
        A(k, [i j]) = 1;
        pairVar(k, :) = [i, j, V(k)];
    end
end
v = A \ V;
clipped = v < 0;
out = struct('sigma', sqrt(max(v, 0))', 'variance', v', ...
    'pairVar', pairVar, 'clipped', clipped');
end
