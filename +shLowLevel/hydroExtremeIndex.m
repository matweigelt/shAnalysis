function [Z, info] = hydroExtremeIndex(X, epochs, opts)
%HYDROEXTREMEINDEX Flood/drought indices per grid cell or basin.
%
%   [Z, INFO] = shLowLevel.hydroExtremeIndex(X, EPOCHS) turns a TWS
%   anomaly stack (grid cells, basin series, or a single series) into
%   a standardized hydrological-extreme index. Three published modes:
%
%   "DSI" (default) - GRACE Drought Severity Index (Zhao et al. 2017):
%       Z_{i,j} = (TWSA_{i,j} - mean_j) / sigma_j per CALENDAR month j
%       - dimensionless, detects drought AND abnormally wet states;
%       validated in the literature against USDM, PDSI, SPEI, NDVI,
%       soil moisture and in-situ groundwater. INFO.category returns
%       the 11-class USDM-style labels (D4..W4).
%   "WSDI" - standardized water storage deficit (Sinha et al. 2017):
%       the climatology-deviation series standardized as one
%       population (not per month).
%   "StorageDeficit" - the Reager & Famiglietti (2009) flood
%       PREDISPOSITION: Sdef(t) = max(X(1..t-1)) - X(t-1), the water
%       the store can still take before reaching its historical
%       maximum - CAUSAL (uses only the past), >= 0, exactly 0 at a
%       running maximum. Small Sdef = saturated basin = flood-prone.
%       With PrecipGrid= the full Flood Potential Index is formed:
%       FPA = P - Sdef, FPI = FPA / max(FPA) per cell (Reager &
%       Famiglietti 2009); without precipitation this function
%       deliberately does NOT call the result FPI.
%
%   Detrend policy, quantified by the Python pre-validation: with a
%   0.5 cm/yr trend (GIA, anthropogenic depletion) a planted
%   exceptional-drought month weakens from -2.41 to -1.44 - OUT of
%   its category - and 45% of the late decade turns spuriously "wet".
%   Detrend = "linear" is therefore the default for DSI/WSDI. For
%   StorageDeficit the Reager original uses NO detrending (the
%   capacity notion is physical); default "none" there, overridable.
%
%   Inputs
%     X       (nLat x nLon x T | K x T | T x 1) double  TWS anomalies
%             [any unit]; grids as from twsChain, basins from
%             shSeries.basinAverage
%     epochs  (T x 1) double  decimal years
%
%   Options
%     Mode ("DSI")        "DSI" | "WSDI" | "StorageDeficit"
%     Detrend ("default") "linear" | "none" | "default"; "default"
%                         resolves to "linear" for DSI/WSDI and to
%                         "none" for StorageDeficit (the physical
%                         Reager convention)
%     Sigma ("std")       "std" | "robust" (1.4826 * MAD - one bad
%                         month corrupts std by 6x, MAD by nothing:
%                         measured in the pre-validation)
%     MinYears (5)        (1 x 1) calendar months sampled by fewer
%                         years get NaN and a warning, never a
%                         quietly unstable sigma
%     PrecipGrid ([])     same shape as X [same unit]: enables the
%                         full FPI in StorageDeficit mode
%
%   Outputs
%     Z    same shape as X  the index (DSI/WSDI z-scores; Sdef in the
%          unit of X; FPI dimensionless in [-inf, 1] when PrecipGrid)
%     info (1 x 1) struct  mode, detrended, clim (12 x cells, DSI),
%          sigma, trendPerYr, category (same shape as Z, int8 -5..+5,
%          DSI only; edges [-2 -1.6 -1.3 -0.8 -0.5 .5 .8 1.3 1.6 2]),
%          categoryNames (11 x 1 string, "D4".."W4"), nPerMonth
%
%   Example
%     [tws, rep] = shLowLevel.twsChain(ser, kn = kn);       % fields
%     [Z, inf1] = shLowLevel.hydroExtremeIndex(tws.grid, tws.epochs);
%     imagesc(inf1.category(:, :, end))    % current drought map
%
%   Error identifiers
%     shLowLevel:hydroExtremeIndex:badSize     epochs vs X mismatch
%     shLowLevel:hydroExtremeIndex:badPrecip   PrecipGrid shape/mode
%
%   References: Zhao, A, Velicogna, Kimball (2017) J Hydrometeorol 18;
%   Reager & Famiglietti (2009) GRL 36, L23402; Reager, Thomas,
%   Famiglietti (2014) Nature Geosci 7; Sinha et al. (2017).
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-13 (v3.14.0).
arguments
    X double
    epochs (:,1) double
    opts.Mode (1,1) string ...
        {mustBeMember(opts.Mode, ["DSI", "WSDI", "StorageDeficit"])} = "DSI"
    opts.Detrend (1,1) string ...
        {mustBeMember(opts.Detrend, ["linear", "none", "default"])} = "default"
    opts.Sigma (1,1) string ...
        {mustBeMember(opts.Sigma, ["std", "robust"])} = "std"
    opts.MinYears (1,1) double {mustBePositive} = 5
    opts.PrecipGrid double = []
end
T = numel(epochs);
sz = size(X);
if sz(end) ~= T && ~(isvector(X) && numel(X) == T)
    error('shLowLevel:hydroExtremeIndex:badSize', ...
        'last dimension of X (%d) must match epochs (%d).', sz(end), T);
end
if isvector(X), X = reshape(X, 1, T); sz = size(X); end
Q = prod(sz(1:end-1));
Xm = reshape(X, Q, T);                       % cells x time
if opts.Detrend == "default"
    if opts.Mode == "StorageDeficit", opts.Detrend = "none";
    else, opts.Detrend = "linear"; end
end
% ---- detrend
trendPerYr = zeros(Q, 1);
Xa = Xm;
if opts.Detrend == "linear"
    A = [ones(T, 1), epochs - mean(epochs)];
    co = A \ Xm';                            % 2 x Q
    Xa = Xm - (A * co)';
    trendPerYr = co(2, :)';
end
mo = max(1, min(12, 1 + floor(mod(epochs, 1) * 12)));
switch opts.Mode
    case {"DSI", "WSDI"}
        clim = nan(Q, 12); sig = nan(Q, 12); nPerMonth = zeros(12, 1);
        Za = nan(Q, T);
        for j = 1:12
            s = mo == j;
            nPerMonth(j) = nnz(s);
            if nnz(s) < opts.MinYears, continue, end
            V = Xa(:, s);
            clim(:, j) = mean(V, 2);
            if opts.Sigma == "robust"
                sig(:, j) = 1.4826 * median(abs(V - median(V, 2)), 2);
            else
                sig(:, j) = std(V, 0, 2);
            end
            Za(:, s) = (V - clim(:, j)) ./ sig(:, j);
        end
        if any(nPerMonth < opts.MinYears)
            warning('shLowLevel:hydroExtremeIndex:shortMonths', ...
                '%d calendar month(s) sampled by fewer than %d years - NaN there.', ...
                nnz(nPerMonth < opts.MinYears), opts.MinYears);
        end
        if opts.Mode == "WSDI"
            % Sinha: standardize the deficit series as one population
            D = Xa - clim(:, mo);
            Za = (D - mean(D, 2, 'omitnan')) ./ std(D, 0, 2, 'omitnan');
        end
        Z = reshape(Za, sz);
        edges = [-inf -2 -1.6 -1.3 -0.8 -0.5 0.5 0.8 1.3 1.6 2 inf];
        cat = int8(discretize(Za, edges)) - 6;         % -5..+5
        info = struct('mode', opts.Mode, 'detrended', opts.Detrend, ...
            'clim', clim, 'sigma', sig, 'trendPerYr', trendPerYr, ...
            'category', reshape(cat, sz), ...
            'categoryNames', ["D4";"D3";"D2";"D1";"D0";"normal"; ...
                "W0";"W1";"W2";"W3";"W4"], ...
            'nPerMonth', nPerMonth);
    case "StorageDeficit"
        runMax = cummax(Xa, 2);
        Sdef = [nan(Q, 1), runMax(:, 1:end-1) - Xa(:, 1:end-1)];
        if ~isempty(opts.PrecipGrid)
            P = reshape(opts.PrecipGrid, [], T);
            if size(P, 1) ~= Q
                error('shLowLevel:hydroExtremeIndex:badPrecip', ...
                    'PrecipGrid must match X in shape.');
            end
            FPA = P - Sdef;
            Z = reshape(FPA ./ max(FPA, [], 2), sz);   % Reager 2009
            what = "FPI";
        else
            Z = reshape(Sdef, sz);
            what = "StorageDeficit";
        end
        info = struct('mode', what, 'detrended', opts.Detrend, ...
            'trendPerYr', trendPerYr, 'note', ...
            "Sdef is causal (past-only), >= 0, exactly 0 at a running maximum");
end
end
