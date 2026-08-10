function P = legendreCached(nmax, latRad)
%LEGENDRECACHED Cached fully normalized Legendre functions.
%
%   P = shx.legendreCached(NMAX, LATRAD) returns shx.legendreALF(NMAX,
%   LATRAD), serving repeated requests for the same (NMAX, LATRAD) from a
%   process-wide cache. Cache hits are verified with an EXACT isequal
%   comparison of the stored latitude vector -- a mismatched latVec can
%   never silently return the wrong P (this removes the documented
%   footgun of shSynthesis' raw 'P' passthrough).
%
%   The cache holds at most 4 entries (FIFO); at nmax=180 with 181
%   latitudes one entry is ~47 MB. Call shx.legendreCached('clear') to
%   empty it.
%
%   Inputs
%     nmax    (1,1) double   maximum degree, or the char 'clear'
%     latRad  (1,:)/(:,1) double  geocentric latitudes [rad]
%   Outputs
%     P       (nmax+1)x(nmax+1)xnumel(latRad) double
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

persistent keys lats vals

if nargin == 1 && (ischar(nmax) || isstring(nmax)) && strcmpi(nmax, 'clear')
    keys = []; lats = {}; vals = {};
    return
end

latRad = latRad(:)';
key = [nmax, numel(latRad), sum(latRad), latRad(1), latRad(end)];

if ~isempty(keys)
    for k = 1:size(keys, 1)
        if isequal(keys(k, :), key) && isequal(lats{k}, latRad)
            P = vals{k};
            return
        end
    end
end

% symmetric-latitude parity trick (exact): recursion only for unique |lat|
[uab, ~, iu] = unique(abs(latRad));
if numel(uab) < numel(latRad)
    Pu = shx.legendreALF(nmax, uab);
    P = Pu(:, :, iu);
    parity = (-1).^((0:nmax)' + (0:nmax));
    P(:, :, latRad < 0) = P(:, :, latRad < 0) .* parity;
else
    P = shx.legendreALF(nmax, latRad);
end

if size(keys, 1) >= 4                      % FIFO eviction
    keys(1, :) = []; lats(1) = []; vals(1) = [];
end
keys(end+1, :) = key;
lats{end+1} = latRad;
vals{end+1} = P;
end
