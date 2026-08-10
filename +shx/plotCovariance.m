function h = plotCovariance(M, idx, opts)
%PLOTCOVARIANCE Covariance/correlation image with order-block boundaries.
%
%   H = shx.plotCovariance(M, IDX) shows a P x P covariance (e.g. from
%   shx.readSINEX with Index=IDX) as an image in shIndex ordering, with
%   the (order, C/S) block boundaries drawn - the one-look check of how
%   good the block-diagonal approximation (tvANS block path, DDK
%   structure) is for THIS matrix.
%
%   Inputs
%     M    (P,P) double     covariance in IDX ordering
%     idx  struct           shx.shIndex
%   Options
%     Type ("correlation")  "correlation" (normalized, [-1,1], divergent
%                           scale) | "covariance" (log10|value|)
%     BlockLines (true)     draw order/cs boundaries
%     ax ([])
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Outputs
%     h          (1 x 1) graphics handle   image of the covariance/correlation structure
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    M double
    idx (1,1) struct
    opts.Type (1,1) string ...
        {mustBeMember(opts.Type, ["correlation","covariance"])} = "correlation"
    opts.BlockLines (1,1) logical = true
    opts.ax = []
end
if isempty(opts.ax), ax = gca; else, ax = opts.ax; end
assert(isequal(size(M), [idx.P idx.P]), 'shx:plotCovariance:badSize', ...
    'M must be P x P for the supplied idx (P = %d).', idx.P);
if opts.Type == "correlation"
    d = sqrt(diag(M)); d(d == 0) = 1;
    img = M ./ (d * d');
    imagesc(ax, img);
    clim(ax, [-1 1]);
    cb = colorbar(ax); cb.Label.String = 'correlation';
else
    img = log10(abs(M)); img(~isfinite(img)) = NaN;
    imagesc(ax, img);
    cb = colorbar(ax); cb.Label.String = 'log_{10}|cov|';
end
axis(ax, 'image');
xlabel(ax, 'parameter index'); ylabel(ax, 'parameter index');
title(ax, sprintf('%s (shIndex order, Lmax=%d)', opts.Type, idx.Lmax));
if opts.BlockLines
    % boundaries where (m, cs) changes
    key = idx.m(:) * 2 + idx.cs(:);
    bnd = find(diff(key) ~= 0) + 0.5;
    hold(ax, 'on');
    for b = bnd(:)'
        xline(ax, b, 'k-', 'Alpha', 0.25);
        yline(ax, b, 'k-', 'Alpha', 0.25);
    end
    hold(ax, 'off');
end
h = ax;
end
