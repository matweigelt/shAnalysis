function h = plotSHCoeffTriangle(C, S, varargin)
%PLOTSHCOEFFTRIANGLE Coefficient-triangle diagnostic plot.
%
%   PLOTSHCOEFFTRIANGLE(C, S) shows |Cnm| (left half) and |Snm| (right
%   half) as a single degree/order "butterfly" triangle image on a log
%   color scale -- the standard first diagnostic for spotting striping,
%   leakage, or bad coefficients before spectral/spatial analysis.
%   x-axis: order m (mirrored: Snm on the LEFT, Cnm on the RIGHT),
%   y-axis: degree n with the apex (n=0) at the TOP - the triangle's
%   peak points upwards. The data aspect ratio is fixed to 1:1, so the
%   axes box is 2:1 (x spans 2*nmax+1 units, y spans nmax+1) and every
%   imagesc cell renders square (v2.3.1).
%
%   Name/value options:
%     'nmax'    max degree to display, default size(C,1)-1
%     'clim'    [min max] log10(|value|) color limits; default auto
%     'ax'      axes handle to plot into (default: current axes)
%     'RefC', 'RefS'  reference coefficients (v2.4): plots the SIGNED
%               difference C-RefC / S-RefS on a linear divergent scale,
%               symmetric about zero - the before/after diagnostic for
%               filter comparisons (e.g. tvANS vs DDK3).
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shx: Claude (Fable 5), 2026-08-07.
%   Outputs
%     h          (1 x 1) graphics handle   butterfly triangle image (abs log10 or signed diff with RefC/RefS)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

p = inputParser;
addParameter(p, 'nmax', size(C,1)-1);
addParameter(p, 'clim', []);
addParameter(p, 'ax', []);
addParameter(p, 'RefC', []);
addParameter(p, 'RefS', []);
parse(p, varargin{:});
nmax = p.Results.nmax;
diffMode = ~isempty(p.Results.RefC);
if diffMode
    assert(~isempty(p.Results.RefS), 'shx:plotSHCoeffTriangle:needRefS', ...
        'RefC requires RefS.');
    C = C - p.Results.RefC(1:size(C,1), 1:size(C,2));
    S = S - p.Results.RefS(1:size(S,1), 1:size(S,2));
end

if isempty(p.Results.ax)
    % MATLAB convention (like plot/imagesc): render into the CURRENT axes,
    % creating a figure only if none exists. Unconditional 'figure' broke
    % subplot workflows: demo figure 2 prepared subplot(1,2,1), triangle
    % opened a fresh window, and the prepared axes stayed empty (v2.2.2).
    ax = gca;
else
    ax = p.Results.ax;
end

if diffMode
    Cv = C(1:nmax+1, 1:nmax+1);
    Sv = S(1:nmax+1, 1:nmax+1);
    tl = tril(true(nmax+1)); Cv(~tl) = NaN;
    tl1 = tl; tl1(:, 1) = false;    % S: n >= m >= 1 INCLUDING sectorals
    Sv(~tl1) = NaN;                 % (tril(...,-1) wrongly blanked S_nn
                                    %  and made the S wing one column
                                    %  narrower in every diff triangle;
                                    %  fixed v2.4.2)
    img = nan(nmax+1, 2*nmax+1);
    img(:, nmax+1:end) = Cv;
    img(:, 1:nmax) = fliplr(Sv(:, 2:end));
else
    Cabs = abs(C(1:nmax+1, 1:nmax+1));
    Sabs = abs(S(1:nmax+1, 1:nmax+1));
    Cabs(Cabs==0) = NaN;
    Sabs(Sabs==0) = NaN;
    % mirrored image: columns -nmax..-1 = Snm (reversed), 0..nmax = Cnm
    img = nan(nmax+1, 2*nmax+1);
    img(:, nmax+1:end) = log10(Cabs);          % m = 0..nmax -> center..right
    img(:, 1:nmax) = fliplr(log10(Sabs(:,2:end))); % m = 1..nmax -> left
end

imagesc(ax, -nmax:nmax, 0:nmax, img);
set(ax, 'YDir', 'reverse');                 % apex n=0 on top: peak upwards
daspect(ax, [1 1 1]);                       % square cells, 2:1 axes box
axis(ax, 'tight');
if diffMode
    half = 128; t = linspace(0, 1, half)';
    bl = [0.10 0.20 0.65]; rd = [0.70 0.05 0.10]; w = [1 1 1];
    colormap(ax, [bl + t.*(w - bl); flipud(rd + t.*(w - rd))]);
    a = max(abs(img(:)), [], 'omitnan');
    if isempty(a) || a == 0, a = 1; end
    clim(ax, [-a, a]);
    cb = colorbar(ax);
    cb.Label.String = 'coefficient difference';
else
    colormap(ax, 'parula');
    cb = colorbar(ax);
    cb.Label.String = 'log_{10}|coefficient|';
end
if ~isempty(p.Results.clim)
    caxis(ax, p.Results.clim); %#ok<CAXIS>
end
xlabel(ax, '\leftarrow order m (S_{nm})      order m (C_{nm}) \rightarrow');
ylabel(ax, 'degree n');
if diffMode
    title(ax, 'Stokes coefficient triangle (signed difference)');
else
    title(ax, 'Stokes coefficient triangle (log_{10}|value|)');
end
h = ax;

end
