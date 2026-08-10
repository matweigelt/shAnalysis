function x = vecFromCS(cnm, snm, idx)
%VECFROMCS Pack lower-triangular Cnm/Snm matrices into a shLowLevel vector.
%   x = shLowLevel.vecFromCS(cnm, snm, idx) with cnm, snm of size
%   (Lmax+1)x(Lmax+1), cnm(n+1,m+1) = Cnm, snm(n+1,m+1) = Snm.
%   Outputs
%     x          (P x 1) double   coefficients gathered in idx ordering
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    cnm double
    snm double
    idx (1,1) struct
end

x  = zeros(idx.P, 1);
iC = idx.cs == 0;
x(iC)  = cnm(sub2ind(size(cnm), idx.n(iC)+1,  idx.m(iC)+1));
x(~iC) = snm(sub2ind(size(snm), idx.n(~iC)+1, idx.m(~iC)+1));
end
