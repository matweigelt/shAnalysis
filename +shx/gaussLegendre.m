function [x, w] = gaussLegendre(n)
%GAUSSLEGENDRE Gauss-Legendre nodes and weights on [-1, 1].
%
%   [X, W] = shx.gaussLegendre(N) via Golub-Welsch (symmetric eigenvalue
%   problem of the Jacobi matrix). Exact for polynomials up to degree
%   2N-1: integral f(t) dt = sum(W .* f(X)).
%
%   Inputs
%     n   (1,1) double   number of nodes
%   Outputs
%     x   (n,1) double   nodes, ascending
%     w   (n,1) double   weights, sum(w) = 2
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    n (1,1) double {mustBeInteger, mustBePositive}
end

k = 1:n-1;
beta = k ./ sqrt(4*k.^2 - 1);
J = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(J);
[x, ord] = sort(diag(D));
w = 2 * V(1, ord)'.^2;
w = w(:);
end
