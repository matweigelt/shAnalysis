function [Cf, Sf] = shDestripe(C, S, varargin)
%SHDESTRIPE Remove correlated (north-south striping) noise from Stokes
%   coefficients, following Swenson & Wahr (2006).
%
%   [CF, SF] = SHDESTRIPE(C, S) removes the smoothly-varying-with-degree
%   component of C(n,m), S(n,m) at each fixed order m (separately for
%   even and odd n, since correlated GRACE errors alias into that
%   parity split), for all orders m >= 'minOrder'. Coefficients at
%   m < minOrder are passed through unchanged.
%
%   Name/value options:
%     'minOrder'     lowest order to filter, default 6 (below this, real
%                    signal dominates and destriping mostly removes signal;
%                    value is somewhat model/application dependent -- verify
%                    against your own signal-to-noise spectrum, cf.
%                    plotSHSpectrum's degree-amplitude vs. error-amplitude)
%     'polyOrder'    polynomial order fit and removed per order/parity,
%                    default 3 (classic Swenson & Wahr cubic fit)
%     'windowLength' if empty (default), fit a single polynomial over the
%                    full available degree range per order/parity (the
%                    original Swenson & Wahr formulation). If set to an
%                    odd integer >= polyOrder+2, use a centered moving
%                    window instead (Duan et al. 2009 / "P3M6"-style
%                    windowed variant, e.g. windowLength=6 with
%                    polyOrder=3 is the commonly used "P3M6" filter).
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shx: Claude (Fable 5), 2026-08-07.
%   Outputs
%     Cf         (nmax+1 x nmax+1) double   destriped cosine coefficients
%     Sf         (nmax+1 x nmax+1) double   destriped sine coefficients
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

p = inputParser;
addParameter(p, 'minOrder', 6);
addParameter(p, 'polyOrder', 3);
addParameter(p, 'windowLength', []);
parse(p, varargin{:});
minOrder = p.Results.minOrder;
polyOrder = p.Results.polyOrder;
windowLength = p.Results.windowLength;

nmax = size(C,1) - 1;
Cf = C;
Sf = S;

for m = minOrder:nmax
    n = (m:nmax)';
    idxEven = n(mod(n,2)==0);
    idxOdd  = n(mod(n,2)==1);

    Cf = filterParity(Cf, idxEven, m, polyOrder, windowLength);
    Cf = filterParity(Cf, idxOdd,  m, polyOrder, windowLength);
    if m > 0
        Sf = filterParity(Sf, idxEven, m, polyOrder, windowLength);
        Sf = filterParity(Sf, idxOdd,  m, polyOrder, windowLength);
    end
end

end

% ------------------------------------------------------------------
function Mat = filterParity(Mat, nIdx, m, polyOrder, windowLength)
%FILTERPARITY Remove polynomial-in-degree trend for one order/parity subset.
if numel(nIdx) <= polyOrder + 1
    return; % not enough points to fit meaningfully, leave unchanged
end
y = Mat(nIdx+1, m+1);

if isempty(windowLength)
    coeff = polyfit(nIdx, y, polyOrder);
    trend = polyval(coeff, nIdx);
else
    trend = localMovingPolyfit(nIdx, y, polyOrder, windowLength);
end

Mat(nIdx+1, m+1) = y - trend;
end

% ------------------------------------------------------------------
function trend = localMovingPolyfit(x, y, polyOrder, windowLength)
%LOCALMOVINGPOLYFIT Centered moving-window polynomial fit, shrinking window
%   near the edges of the series (still >= polyOrder+1 points required).
N = numel(x);
half = floor(windowLength/2);
trend = zeros(N,1);
for k = 1:N
    lo = max(1, k-half);
    hi = min(N, k+half);
    while (hi-lo+1) < polyOrder+1 && (lo > 1 || hi < N)
        if lo > 1, lo = lo - 1; end
        if hi < N, hi = hi + 1; end
    end
    coeff = polyfit(x(lo:hi), y(lo:hi), polyOrder);
    trend(k) = polyval(coeff, x(k));
end
end
