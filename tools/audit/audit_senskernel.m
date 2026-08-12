function audit_senskernel
%AUDIT_SENSKERNEL Reproduce the 2-16% claim: tailored vs Gaussian at
%   matched noise & gain=1, across FarField weightings. Audit script,
%   Claude (Fable 5), 2026-08-12.
idx = shLowLevel.shIndex(30, MinDegree = 0);
cap = @(la, lo) double(acosd(min(1, max(-1, sind(15)*sind(la) + ...
    cosd(15)*cosd(la).*cosd(lo - 300)))) <= 10);
% default degree-dependent GRACE-like noise (function's own default)
[~, i0] = shLowLevel.sensitivityKernel(idx, cap, Alpha = 0);
kEx = i0.kExact; N = defaultNoiseVec(idx);        % rebuild same default
weights = {"1/(n+1) (default)", 1./(idx.n+1);
           "flat",              ones(idx.P,1);
           "high-degree n",     max(idx.n,1)};
fprintf('%-20s %-9s %-9s %-9s %-7s\n','FarField','noiseTgt','leakG','leakT','margin');
for wsp = 1:size(weights,1)
    M = weights{wsp,2};
    for rGauss = [400 600 800]
        % Gaussian competitor: smoothed indicator, renormalised to gain 1
        w = shLowLevel.shGaussianWeights(idx.Lmax, rGauss); w = w(:);
        kG = kEx .* w(idx.n+1);
        kG = kG * (kEx'*kEx) / (kG'*kEx);          % gain = 1
        noiseTgt = sqrt(sum(N .* kG.^2));
        leakG = sqrt(sum(M .* (kG - kEx).^2));
        % tailored at MATCHED noise: bisect Alpha
        aL = 1e-6; aH = 1e8;
        for it = 1:60
            a = sqrt(aL*aH);
            [kT, iT] = shLowLevel.sensitivityKernel(idx, cap, Alpha=a, ...
                FarField=M, Noise=sqrt(N));   % vector input = SIGMAS
            if iT.noise > noiseTgt, aL = a; else, aH = a; end
        end
        leakT = sqrt(sum(M .* (kT - kEx).^2));
        assert(abs(iT.gain - 1) < 1e-8);
        fprintf('%-20s %-9.3g %-9.3g %-9.3g %+6.1f%%  (r=%d, alpha=%.3g, noiseT=%.3g)\n', ...
            weights{wsp,1}, noiseTgt, leakG, leakT, 100*(leakT-leakG)/leakG, ...
            rGauss, a, iT.noise);
    end
end
end

function N = defaultNoiseVec(idx)
% mirror of the function's default (sensitivityKernel.m:124): variance
% vector with GRACE-like growth; shape only, Alpha absorbs scale
N = (1 + (idx.n(:) / 8).^3).^2;
end
