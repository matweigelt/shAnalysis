function h = plotSHMap(grid, latDeg, lonDeg, opts)
%PLOTSHMAP Global map plot with coastlines and a sane divergent scale.
%
%   H = shLowLevel.plotSHMap(GRID, LAT, LON) renders a synthesized field the way
%   a GRACE map should look by default: divergent blue-white-red colormap
%   symmetric about zero, robust color limits (98th percentile of |value|),
%   coastline overlay (base MATLAB 'coastlines'), correct aspect. The
%   convenience method g.map(...) synthesizes and plots in one call.
%
%   Inputs
%     grid    (nlat,nlon) double
%     latDeg  (1,nlat) double [deg]
%     lonDeg  (1,nlon) double [deg]
%   Options
%     Projection ("plate")  "plate" (equirectangular) | "hammer"
%                           (equal-area, analytic - no Mapping Toolbox)
%     Coast (true)          coastline overlay
%     CLim ([])             color limits; default symmetric robust
%     Colormap ("divergent") "divergent" (blue-white-red, built
%                           programmatically) | any colormap name
%     Title (""), Units ("")  annotations
%     ax ([])               target axes (default: current axes)
%   Outputs
%     h  (1,1) graphics handle  axes of the map plot
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    grid double
    latDeg (1,:) double
    lonDeg (1,:) double
    opts.Projection (1,1) string ...
        {mustBeMember(opts.Projection, ["plate","hammer"])} = "plate"
    opts.Coast (1,1) logical = true
    opts.CLim double = []
    opts.Colormap (1,1) string = "divergent"
    opts.Title (1,1) string = ""
    opts.Units (1,1) string = ""
    opts.ax = []
end
if isempty(opts.ax), ax = gca; else, ax = opts.ax; end
assert(isequal(size(grid), [numel(latDeg), numel(lonDeg)]), ...
    'shLowLevel:plotSHMap:badSize', 'GRID must be nlat x nlon matching LAT/LON.');

cl = opts.CLim;
if isempty(cl)
    a = shLowLevel.pctile(abs(grid), 98);
    if a == 0, a = 1; end
    cl = [-a, a];
end

% wrap longitudes to a monotonic axis centered on the data
lonP = lonDeg(:)';
if opts.Projection == "plate"
    imagesc(ax, lonP, latDeg, grid);
    set(ax, 'YDir', 'normal');
    daspect(ax, [1 1 1]);
    xlabel(ax, 'longitude [deg]'); ylabel(ax, 'latitude [deg]');
else
    % Hammer equal-area: analytic forward projection
    [LON, LAT] = meshgrid(lonP, latDeg);
    [xh, yh] = hammer(LAT, LON);
    surface(ax, xh, yh, zeros(size(xh)), grid, 'EdgeColor', 'none');
    view(ax, 2);
    axis(ax, 'equal', 'off');
end
if ~isempty(cl), clim(ax, cl); end
if opts.Colormap == "divergent"
    colormap(ax, divergentMap(256));
else
    colormap(ax, char(opts.Colormap));
end
cb = colorbar(ax);
if strlength(opts.Units) > 0, cb.Label.String = char(opts.Units); end
if strlength(opts.Title) > 0, title(ax, opts.Title); end

if opts.Coast
    try
        cst = load('coastlines');
        cla_ = cst.coastlat; clo_ = mod(cst.coastlon, 360);
        % break wrap-around jumps so no spurious horizontal lines appear
        jmp = [false; abs(diff(clo_)) > 180];
        cla_(jmp) = NaN; clo_(jmp) = NaN;
        hold(ax, 'on');
        if opts.Projection == "plate"
            plot(ax, clo_, cla_, 'k-', 'LineWidth', 0.5);
        else
            [cx, cy] = hammer(cla_, clo_);
            plot(ax, cx, cy, 'k-', 'LineWidth', 0.5);
        end
        hold(ax, 'off');
    catch
        % coastlines.mat unavailable in some deployments: map still valid
    end
end
h = ax;
end

function [x, y] = hammer(latDeg, lonDeg)
%HAMMER Analytic Hammer-Aitoff projection, lon centered at 180 deg.
phi = deg2rad(latDeg);
lam = deg2rad(mod(lonDeg - 180, 360) - 180);       % [-pi, pi] about 180E
z = sqrt(1 + cos(phi) .* cos(lam/2));
x = 2 * sqrt(2) * cos(phi) .* sin(lam/2) ./ z;
y = sqrt(2) * sin(phi) ./ z;
end

function cmap = divergentMap(n)
%DIVERGENTMAP Blue-white-red, programmatic (no toolbox dependency).
half = floor(n/2);
t = linspace(0, 1, half)';
blue = [0.10 0.20 0.65]; red = [0.70 0.05 0.10]; w = [1 1 1];
cmap = [blue + t .* (w - blue); flipud(red + t .* (w - red))];
end
