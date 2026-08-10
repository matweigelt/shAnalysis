function [cnm, snm] = csFromVec(x, idx)
%CSFROMVEC Unpack a shLowLevel vector into lower-triangular Cnm/Snm matrices.
%   [cnm, snm] = shLowLevel.csFromVec(x, idx). Coefficients outside the index
%   range (n < MinDegree) are zero.
%   Outputs
%     C          (nmax+1 x nmax+1) double   cosine coefficients scattered from x
%     S          (nmax+1 x nmax+1) double   sine coefficients scattered from x
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    x (:,1) double
    idx (1,1) struct
end

L = idx.Lmax;
cnm = zeros(L+1); snm = zeros(L+1);
iC = idx.cs == 0;
cnm(sub2ind(size(cnm), idx.n(iC)+1,  idx.m(iC)+1))  = x(iC);
snm(sub2ind(size(snm), idx.n(~iC)+1, idx.m(~iC)+1)) = x(~iC);
end
