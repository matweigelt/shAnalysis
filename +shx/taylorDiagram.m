function h = taylorDiagram(refStd, stds, corrs, opts)
%TAYLORDIAGRAM Taylor diagram (std / correlation / centered RMSD).
%
%   H = shx.taylorDiagram(REFSTD, STDS, CORRS) draws the standard
%   one-figure comparison of solutions against a reference: radius is
%   the (weighted) standard deviation, the azimuth is acos(correlation),
%   and by the law of cosines the distance to the reference point equals
%   the centered RMSD - so pattern amplitude, pattern correlation and
%   pattern error are read off a single plot. Feed it the .stdA/.stdB/
%   .corr fields of shx.spatialStats, or temporal statistics from
%   shx.compareSeries. Base MATLAB graphics only.
%
%   Options
%     Labels ([])         string array, one per point
%     Title ("Taylor diagram")
%     Normalize (false)   divide all stds by REFSTD (radius 1 = reference)
%
%   Outputs
%     h          (1,1) graphics handle  the axes handle
%
%   Example
%     st = shx.spatialStats(grid, 0.8 * grid, lat, lon);
%     shx.taylorDiagram(st.stdA, st.stdB, st.corr, Labels = "damped");
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.6.0).
arguments
    refStd (1,1) double {mustBePositive}
    stds (1,:) double
    corrs (1,:) double
    opts.Labels string = strings(1, 0)
    opts.Title (1,1) string = "Taylor diagram"
    opts.Normalize (1,1) logical = false
end
if numel(stds) ~= numel(corrs)
    error('shx:taylorDiagram:sizeMismatch', ...
        'STDS and CORRS must have the same length.');
end
s = stds; r0 = refStd;
if opts.Normalize, s = s / refStd; r0 = 1; end
rmax = 1.35 * max([r0, s(:)']);
h = gca; hold(h, 'on'); axis(h, 'equal');
% std arcs
th = linspace(0, pi/2, 90);
for rr = linspace(rmax/4, rmax, 4)
    plot(h, rr*cos(th), rr*sin(th), ':', 'Color', [0.75 0.75 0.75]);
end
% correlation rays + labels
for rc = [0.2 0.4 0.6 0.7 0.8 0.9 0.95 0.99]
    a = acos(rc);
    plot(h, [0 rmax*cos(a)], [0 rmax*sin(a)], ':', 'Color', [0.85 0.85 0.85]);
    text(h, 1.03*rmax*cos(a), 1.03*rmax*sin(a), sprintf('%.2g', rc), ...
        'FontSize', 8, 'HorizontalAlignment', 'left');
end
% centered-RMSD arcs about the reference point
for rr = linspace(r0/2, 1.5*r0, 3)
    xx = r0 + rr*cos(th + pi/2); yy = rr*sin(th + pi/2);
    keep = hypot(xx, yy) <= rmax & xx >= 0;
    plot(h, xx(keep), yy(keep), '--', 'Color', [0.85 0.7 0.7]);
end
plot(h, r0, 0, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
text(h, r0, -0.06*rmax, 'ref', 'HorizontalAlignment', 'center');
cols = lines(numel(s));
for k = 1:numel(s)
    a = acos(min(max(corrs(k), -1), 1));
    plot(h, s(k)*cos(a), s(k)*sin(a), 'o', 'Color', cols(k, :), ...
        'MarkerFaceColor', cols(k, :), 'MarkerSize', 7);
    if k <= numel(opts.Labels)
        text(h, s(k)*cos(a) + 0.02*rmax, s(k)*sin(a), opts.Labels(k), ...
            'FontSize', 8);
    end
end
xlim(h, [0 1.12*rmax]); ylim(h, [0 1.12*rmax]);
xlabel(h, 'standard deviation'); title(h, opts.Title);
text(h, 0.72*rmax, 1.05*rmax, 'correlation', 'FontSize', 8);
hold(h, 'off');
end
