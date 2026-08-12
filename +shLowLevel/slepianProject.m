function [a, rec] = slepianProject(Cs, Ss, G, idx)
%SLEPIANPROJECT Project coefficient sets onto a Slepian basis (and back).
%
%   [A, REC] = shLowLevel.slepianProject(CS, SS, G, IDX) projects one or
%   many coefficient sets onto the Slepian tapers G from
%   shLowLevel.slepianBasis: A = G' * X where X stacks the (C, S)
%   coefficients in IDX ordering. This is the missing application half
%   of the concentration problem (roadmap item 8): a regional time
%   series analysis estimates only the ~Shannon leading Slepian
%   coefficients instead of all P spherical-harmonic ones - the
%   well-posed alternative to Kaula-regularized regional least squares.
%   REC returns the reconstruction G*A repacked to (C, S) arrays for
%   direct synthesis.
%
%   Inputs
%     Cs   (n+1 x n+1 x T) cosine coefficients (T = 1 for a single set)
%     Ss   (n+1 x n+1 x T) sine coefficients
%     G    (P x K) Slepian taper matrix from slepianBasis (any K)
%     idx  (1 x 1) struct  the SAME shLowLevel.shIndex the basis was
%          built with (P must match size(G, 1))
%
%   Outputs
%     a    (K x T) Slepian-domain coefficient time series
%     rec  (1 x 1) struct with fields
%       Cs (n+1 x n+1 x T) reconstructed cosine coefficients
%       Ss (n+1 x n+1 x T) reconstructed sine coefficients
%
%   Example
%     idx = shLowLevel.shIndex(30, MinDegree = 0);
%     [G, lam, info] = shLowLevel.slepianBasis(idx, @(la, lo) la > 60);
%     K = round(info.shannon);
%     a = shLowLevel.slepianProject(ts.Cs, ts.Ss, G(:, 1:K), idx);
%     plot(ts.epochs, a')          % regional modes vs time
%
%   Error identifiers
%     shLowLevel:slepianProject:sizeMismatch  G rows do not match idx.P
%       or the coefficient arrays are smaller than idx.Lmax
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.9.0).
arguments
    Cs double
    Ss double
    G (:,:) double
    idx (1,1) struct
end
P = idx.P;
if size(G, 1) ~= P
    error('shLowLevel:slepianProject:sizeMismatch', ...
        'G has %d rows but idx.P = %d.', size(G, 1), P);
end
if size(Cs, 1) < idx.Lmax + 1
    error('shLowLevel:slepianProject:sizeMismatch', ...
        'coefficients reach degree %d but idx.Lmax = %d.', ...
        size(Cs, 1) - 1, idx.Lmax);
end
T = size(Cs, 3);
li = sub2ind([size(Cs, 1), size(Cs, 2)], idx.n(:) + 1, idx.m(:) + 1);
isC = idx.cs(:) == 0;
X = zeros(P, T);
for t = 1:T
    Ck = Cs(:, :, t); Sk = Ss(:, :, t);
    v = zeros(P, 1);
    v(isC) = Ck(li(isC));
    v(~isC) = Sk(li(~isC));
    X(:, t) = v;
end
a = G' * X;
if nargout > 1
    XR = G * a;
    nmax = size(Cs, 1) - 1;
    rec = struct('Cs', zeros(nmax+1, nmax+1, T), 'Ss', zeros(nmax+1, nmax+1, T));
    for t = 1:T
        Ck = zeros(nmax+1); Sk = zeros(nmax+1);
        Ck(li(isC)) = XR(isC, t);
        Sk(li(~isC)) = XR(~isC, t);
        rec.Cs(:, :, t) = Ck; rec.Ss(:, :, t) = Sk;
    end
end
end
