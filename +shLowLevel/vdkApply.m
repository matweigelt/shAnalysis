function [xf, info] = vdkApply(x, N, degOfPar, sigM, opts)
%VDKAPPLY Apply the VDK/VADER decorrelation filter to one solution.
%
%   XF = shLowLevel.vdkApply(X, N, DEGOFPAR, SIGM) computes
%
%       xf = (N + alpha * M) \ (N * x)      (Horvath et al. 2018, Eq. 1)
%
%   with the monthly normal-equation matrix N (the inverse error
%   covariance, from the ITSG SINEX via shLowLevel.readSINEX with
%   Index = idx), the diagonal inverse-signal-variance matrix
%   M = diag(1 ./ sigM(l)^2) evaluated per parameter degree, and the
%   filter strength alpha. The filter matrix is never formed: one
%   Cholesky factorization of (N + alpha M) and a solve - O(P^3) per
%   month, seconds at n96 (P = 9405).
%
%   Relation to tvANSFilter, algebraically: with S = M^-1 this is
%   identical to the Wiener form S (S + alpha N^-1)^-1 the tvANS chain
%   uses - both are [I + alpha N^-1 M]^-1 (Python-prevalidated to
%   1e-15). The difference is the INPUT: VDK uses the formal monthly
%   N whose STRUCTURE changes month to month (orbit, repeat cycles,
%   instrument state), tvANS estimates one empirical structure from
%   the series and scales it per month (VCE). The tvANS
%   one-eigendecomposition efficiency trick requires that stationary
%   structure and is therefore unavailable here. Where monthly
%   covariances exist, VDK is the stronger filter (15% median
%   cumulative geoid error reduction, outlier months an order of
%   magnitude - Horvath et al. 2018 Table 2); tvANS remains the tool
%   for series without released covariances.
%
%   Inputs
%     x        (P x 1) double  coefficient vector, same ordering as N
%     N        (P x P) double  normal-equation matrix (symmetric
%              positive definite)
%     degOfPar (P x 1) double  spherical-harmonic degree of every
%              parameter (idx.n from shLowLevel.shIndex)
%     sigM     (P x 1 | 1 x 2) double  signal sigma per parameter, OR
%              [a, b] evaluating sigma = a * l^b per parameter (one
%              row of signalVarianceKaula's AB)
%
%   Options
%     Alpha (1)  (1 x 1) filter strength; larger = smoother (paper
%                Table 1 maps alpha to mean Gaussian radii for the
%                ITSG-Grace2014 case - the mapping is series-specific)
%
%   Outputs
%     xf   (P x 1) double  filtered coefficient vector
%     info (1 x 1) struct  alpha, P, cholOK (true when (N + alpha M)
%          factorized cleanly)
%
%   Example
%     idx = shLowLevel.shIndex(96, MinDegree = 2);
%     snx = shLowLevel.readSINEX(f, Index = idx);
%     xf = shLowLevel.vdkApply(snx.x, snx.M, idx.n, ab(mo, :), Alpha = 1);
%
%   Error identifiers
%     shLowLevel:vdkApply:badSize  x, N, degOfPar inconsistent
%     shLowLevel:vdkApply:notSPD   (N + alpha M) not positive definite
%
%   Reference: Horvath, Murboeck, Pail, Horwath (2018), Geosciences 8,
%   323, doi:10.3390/geosciences8090323; Kusche (2007), J Geod 81.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.13.0).
arguments
    x (:,1) double
    N (:,:) double
    degOfPar (:,1) double
    sigM double
    opts.Alpha (1,1) double {mustBeNonnegative} = 1
end
P = numel(x);
if ~isequal(size(N), [P, P]) || numel(degOfPar) ~= P
    error('shLowLevel:vdkApply:badSize', ...
        'x (%d), N (%dx%d) and degOfPar (%d) are inconsistent.', ...
        P, size(N), numel(degOfPar));
end
if numel(sigM) == 2 && isrow(sigM)
    s = sigM(1) * degOfPar.^sigM(2);
else
    s = sigM(:);
    if numel(s) ~= P
        error('shLowLevel:vdkApply:badSize', ...
            'sigM must be [a, b] or P x 1 (got %d for P = %d).', ...
            numel(s), P);
    end
end
mDiag = 1 ./ s.^2;
B = N + opts.Alpha * spdiags(mDiag, 0, P, P);
[R, flag] = chol(B);
if flag ~= 0
    error('shLowLevel:vdkApply:notSPD', ...
        '(N + alpha M) is not positive definite (chol flag %d).', flag);
end
xf = R \ (R' \ (N * x));
info = struct('alpha', opts.Alpha, 'P', P, 'cholOK', true);
end
