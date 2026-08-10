function h = plotBasinSeries(tYears, c, sigma, opts)
%PLOTBASINSERIES Basin time series with 1-sigma band, gaps, trend.
%
%   H = shLowLevel.plotBasinSeries(TYEARS, C, SIGMA) draws the standard basin-
%   average figure: shaded 1-sigma band (fill), data line, the GRACE <->
%   GRACE-FO gap (or any gap > GapThreshold (0.3) years) hatched grey, and an
%   annotated trend fitted with the toolbox climatology design (bias/
%   trend/annual/semi-annual, AR(1)-corrected sigma).
%
%   Inputs
%     tYears (T,1) double   epochs [decimal years]
%     c      (T,1) double   basin averages (e.g. out.c from
%                           shLowLevel.basinDeconvolve)
%     sigma  (T,1) double   1-sigma; [] for no band
%   Options
%     Label ("basin average"), Units (""), GapThreshold (0.3) [yr],
%     Trend (true)  fit + annotate trend +/- sigma (AR(1)-corrected),
%     ax ([])
%   Outputs
%     h  (1,1) graphics handle  axes of the basin series plot
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    tYears (:,1) double
    c (:,1) double
    sigma double = []
    opts.Label (1,1) string = "basin average"
    opts.Units (1,1) string = ""
    opts.GapThreshold (1,1) double = 0.3
    opts.Trend (1,1) logical = true
    opts.ax = []
end
if isempty(opts.ax), ax = gca; else, ax = opts.ax; end
[tYears, ord] = sort(tYears); c = c(ord);
if ~isempty(sigma), sigma = sigma(ord); end
hold(ax, 'on');
yl = [min(c - abs(max(c) - min(c))*0.1), max(c + abs(max(c) - min(c))*0.1)];
if ~isempty(sigma)
    yl = [min(c - 2*sigma), max(c + 2*sigma)];
end
% gaps first (background)
dg = find(diff(tYears) > opts.GapThreshold);
for k = dg(:)'
    fill(ax, [tYears(k), tYears(k+1), tYears(k+1), tYears(k)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
if ~isempty(sigma)
    fill(ax, [tYears; flipud(tYears)], ...
        [c - sigma; flipud(c + sigma)], [0.55 0.70 0.90], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', '1\sigma');
end
plot(ax, tYears, c, '.-', 'Color', [0.05 0.25 0.55], ...
    'DisplayName', char(opts.Label));
if opts.Trend && numel(tYears) >= 8
    [~, ~, coef, ~, cs] = shLowLevel.fitDeterministicModel(c', tYears, ...
        ARCorrect = true);
    tt = linspace(tYears(1), tYears(end), 200)';
    t0 = mean(tYears);
    plot(ax, tt, coef(1) + coef(2)*(tt - t0), 'r-', 'LineWidth', 1, ...
        'DisplayName', sprintf('trend %.3g \\pm %.2g %s/yr', ...
        coef(2), cs(2), char(opts.Units)));
end
hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, 'epoch [yr]');
if strlength(opts.Units) > 0, ylabel(ax, char(opts.Units)); end
legend(ax, 'Location', 'best');
if ~isempty(sigma) && yl(2) > yl(1), ylim(ax, yl); end
h = ax;
end
