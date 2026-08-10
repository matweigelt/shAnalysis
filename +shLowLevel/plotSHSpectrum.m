function h = plotSHSpectrum(specs, varargin)
%PLOTSHSPECTRUM Log-log degree-amplitude spectrum plot.
%
%   PLOTSHSPECTRUM(SPEC) plots SPEC.degree vs SPEC.degAmplitude (from
%   SHDEGREERMS) on a log-log scale, in mm geoid-height-equivalent.
%
%   PLOTSHSPECTRUM({SPEC1, SPEC2, ...}, 'names', {'name1','name2',...})
%   overlays multiple spectra for comparison (e.g. successive GRACE months,
%   or signal vs formal error).
%
%   Name/value options:
%     'names'   cell array of legend labels
%     'units'   'mm' (default) or 'm' for the y-axis
%     'field'   which spec field to plot, default 'degAmplitude'
%               (also accepts 'errAmplitude' to overlay error spectra)
%     'ax'      axes handle to plot into (default: current axes)
%     'Quantity' (v2.4.1) "amplitude" (default, [m]) | "rms" | "variance"
%               | "cumamplitude" | "cumrms" | "cumvariance" - which
%               spectral quantity to draw; error overlays ride along
%               when the spec carries the matching err* field
%     Domain is inferred from spec.domain: pass shLowLevel.shOrderRMS output
%               for order-domain plots (striping axis); mixing domains
%               in one call errors. X-axis is LINEAR since v2.4.1.
%     'Kaula'   scalar K (v2.4): overlay the Kaula-rule reference degree
%               amplitude R*sqrt(2n+1)*K/n^2 (typical K ~ 1e-5) as a
%               dashed grey line - the sanity yardstick for signal decay
%     'MarkCrossover' true/false (v2.4): with exactly TWO spectra, mark
%               the signal/error crossover degree found by
%               shLowLevel.shSpectralCrossover with a vertical line
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shLowLevel: Claude (Fable 5),
%   2026-08-07 (Kaula + crossover overlays v2.4).
%   Outputs
%     h          (1 x 1) graphics handle   log-log spectrum axes
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

if ~iscell(specs)
    specs = {specs};
end

p = inputParser;
addParameter(p, 'names', {});
addParameter(p, 'units', 'mm');
addParameter(p, 'field', 'degAmplitude');
addParameter(p, 'ax', []);
addParameter(p, 'Kaula', []);
addParameter(p, 'MarkCrossover', false);
addParameter(p, 'Quantity', 'amplitude');
parse(p, varargin{:});
names = p.Results.names;
units = p.Results.units;
field = p.Results.field;

if isempty(p.Results.ax)
    % MATLAB convention (like plot/imagesc): render into the CURRENT axes,
    % creating a figure only if none exists. Unconditional 'figure' broke
    % subplot workflows: demo figure 2 prepared subplot(1,2,1), triangle
    % opened a fresh window, and the prepared axes stayed empty (v2.2.2).
    ax = gca;
else
    ax = p.Results.ax;
end

% ---- quantity/domain resolution (v2.4.1)
% back-compat: an explicit legacy 'field' wins over Quantity
fieldExplicit = ~any(strcmp(p.UsingDefaults, 'field'));
qty = lower(char(p.Results.Quantity));
dom = 'degree';
if isfield(specs{1}, 'domain'), dom = specs{1}.domain; end
for k = 2:numel(specs)
    dk = 'degree';
    if isfield(specs{k}, 'domain'), dk = specs{k}.domain; end
    assert(strcmp(dk, dom), 'shLowLevel:plotSHSpectrum:mixedDomains', ...
        'All spectra must share one domain (degree or order).');
end
if strcmp(dom, 'order'), pre = 'ord'; xf = 'order';
else, pre = 'deg'; xf = 'degree'; end
switch qty
    case 'amplitude',    field = [pre 'Amplitude']; ef = 'errAmplitude';
        ylab = sprintf('%s amplitude [%%s]', dom); isAmp = true;
    case 'rms',          field = [pre 'RMS'];       ef = 'errRMS';
        ylab = sprintf('%s RMS [-]', dom); isAmp = false;
    case 'variance',     field = [pre 'Variance'];  ef = 'errVariance';
        ylab = sprintf('%s variance [-]', dom); isAmp = false;
    case 'cumamplitude', field = 'cumAmplitude';    ef = '';
        ylab = sprintf('cumulative %s amplitude [%%s]', dom); isAmp = true;
    case 'cumrms',       field = 'cumRMS';          ef = '';
        ylab = sprintf('cumulative %s RMS [-]', dom); isAmp = false;
    case 'cumvariance',  field = 'cumVariance';     ef = '';
        ylab = sprintf('cumulative %s variance [-]', dom); isAmp = false;
    otherwise
        error('shLowLevel:plotSHSpectrum:badQuantity', 'Unknown Quantity: %s', qty);
end
if fieldExplicit
    field = p.Results.field;                 % legacy override
    ef = ''; isAmp = contains(lower(field), 'amplitude');
    ylab = [field ' [%s]'];
end
scale = 1;
if isAmp && strcmpi(units, 'mm'), scale = 1000; end

% error curves ride along when present and matching (one per spec)
plotSpecs = {}; plotNames = {};
for k = 1:numel(specs)
    plotSpecs{end+1} = specs{k}; %#ok<AGROW>
    if k <= numel(names), plotNames{end+1} = names{k}; else, plotNames{end+1} = ''; end %#ok<AGROW>
end

colors = lines(numel(plotSpecs));
hold(ax, 'on');
h = gobjects(0, 1);
hn = {};
for k = 1:numel(plotSpecs)
    s = plotSpecs{k};
    assert(isfield(s, field), 'shLowLevel:plotSHSpectrum:missingField', ...
        'Spectrum %d has no field %s (wrong Quantity/Domain?).', k, field);
    y = s.(field) * scale;
    n = s.(xf);
    valid = y > 0;            % log y-scale needs y > 0
    if strcmp(dom, 'degree')
        valid = valid & n >= 1;   % degree 0 (total mass) stays off the
    end                           % spectrum, as in v2.x; order m=0 plots
    h(end+1, 1) = plot(ax, n(valid), y(valid), '-o', 'Color', colors(k,:), ...
        'MarkerSize', 3, 'LineWidth', 1.2); %#ok<AGROW>
    if ~isempty(plotNames{k}), hn{end+1} = plotNames{k}; ...
    else, hn{end+1} = sprintf('spectrum %d', k); end %#ok<AGROW>
    if ~isempty(ef) && isfield(s, ef)
        ye = s.(ef) * scale;
        ve = ye > 0;
        h(end+1, 1) = plot(ax, n(ve), ye(ve), '--', 'Color', colors(k,:), ...
            'LineWidth', 1.0); %#ok<AGROW>
        hn{end+1} = [hn{end}, ' (error)']; %#ok<AGROW>
    end
end
names = hn;
set(ax, 'XScale', 'linear', 'YScale', 'log');   % linear x (v2.4.1)
grid(ax, 'on');
xlabel(ax, sprintf('Spherical harmonic %s', xf));
ylabel(ax, sprintf(ylab, units));
if ~isempty(p.Results.Kaula)
    assert(strcmp(dom, 'degree') && isAmp && ~startsWith(field, 'cum'), ...
        'shLowLevel:plotSHSpectrum:kaulaDomain', ...
        'Kaula overlay applies to degree-domain amplitude plots only.');
    s1 = specs{1};
    nk = s1.degree(s1.degree >= 2);
    R = s1.R;
    yk = R .* sqrt(2*nk + 1) .* p.Results.Kaula ./ nk.^2 * scale;
    hk = plot(ax, nk, yk, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
    h(end+1) = hk;
    names{end+1} = sprintf('Kaula K=%.0e', p.Results.Kaula);
end
if p.Results.MarkCrossover
    % shSpectralCrossover takes ONE spec carrying degAmplitude AND
    % errAmplitude (shDegreeRMS with sigmaC/sigmaS); with two separate
    % spectra, fuse them (first = signal, second = error)
    try
        sc = specs{1};
        if numel(specs) == 2 && ~isfield(sc, 'errAmplitude')
            sc.errAmplitude = specs{2}.(field);
        end
        [~, nc] = shLowLevel.shSpectralCrossover(sc);
        if isfinite(nc)
            xline(ax, nc, ':', sprintf('crossover n=%.1f', nc), ...
                'Color', [0.3 0.3 0.3]);
        end
    catch
        % no crossover / incompatible specs: overlay silently omitted
    end
end
if ~isempty(names)
    legend(ax, h, names, 'Location', 'northeast');
end
title(ax, 'Spectral (degree-domain) analysis');
hold(ax, 'off');

end
