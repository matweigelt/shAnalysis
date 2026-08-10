function [Cf, Sf] = applyDDK(C, S, W)
%APPLYDDK Apply an order-block-diagonal (DDK-style) filter to coefficients.
%
%   [CF, SF] = shLowLevel.applyDDK(C, S, W) applies the filter container from
%   shLowLevel.readDDK: within each (order m, cs) block, the degree vector
%   c(n1..n2) is replaced by W.blocks(k).M * c. Degrees present in C but
%   not covered by any block pass through UNCHANGED (typically n < 2).
%   User-facing wrappers: g.applyDDK(W), ts.applyDDK(W).
%
%   Inputs
%     C, S  (n1,n1) double   coefficient matrices, C(n+1,m+1)
%     W     struct           from shLowLevel.readDDK (nmax, blocks(m,cs,n,M))
%   Outputs
%     Cf         (nmax+1 x nmax+1) double   DDK-filtered cosine coefficients
%     Sf         (nmax+1 x nmax+1) double   DDK-filtered sine coefficients
%
%   Filter blocks reaching beyond the field's nmax are TRUNCATED to the
%   available degrees (submatrix of the gain; v2.5.1) - the standard
%   practical application of the released Lmax-120 Wbd filters to
%   n60/n96 GRACE fields. Fields with nmax >= the block degree are
%   filtered exactly as before.
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    C double
    S double
    W (1,1) struct
end
n1 = size(C, 1);
Cf = C; Sf = S;
for k = 1:numel(W.blocks)
    b = W.blocks(k);
    if b.m + 1 > n1, continue; end
    % v2.5.1: filter blocks reaching beyond the field's nmax are
    % TRUNCATED to the available degrees (submatrix of the gain) - the
    % standard practical application of the released Lmax-120 Wbd
    % filters to n60/n96 GRACE fields. For fields with nmax >= the
    % block degree the result is unchanged.
    keep = b.n + 1 <= n1;
    if ~any(keep), continue; end
    ii = b.n(keep) + 1;
    M = b.M(keep, keep);
    if b.cs == 0
        Cf(ii, b.m+1) = M * C(ii, b.m+1);
    else
        Sf(ii, b.m+1) = M * S(ii, b.m+1);
    end
end
end
