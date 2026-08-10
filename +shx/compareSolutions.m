function [rep, h] = compareSolutions(g1, g2, opts)
%COMPARESOLUTIONS Full spectral + spatial comparison of two solutions.
%
%   REP = shx.compareSolutions(G1, G2) compares two shCoefficients on a
%   COMMON BASIS - both are truncated to the smaller nmax and G2 is
%   rescaled to G1's GM/R if they differ (otherwise the comparison
%   measures processing chains, not solutions) - and returns the
%   standard metric set:
%     spectral: difference degree amplitude vs both signal spectra,
%               per-degree correlation, crossover degree of agreement
%               (shx.diffSpectrum)
%     spatial:  area-weighted bias/RMSD/centered RMSD/pattern
%               correlation/std ratio of the synthesized quantity
%               (shx.spatialStats; Taylor-ready)
%     errors:   chi^2/dof of the coefficient differences against the
%               combined formal sigmas, when both carry sigmas - this
%               tests error REALISM, a separate question from agreement
%   [REP, H] = shx.compareSolutions(..., Plot = true) adds a 4-panel
%   figure: spectra, degree correlation, difference map, and the order
%   RMS of the difference (stripe and alias differences live at
%   specific orders).
%
%   Options
%     LatDeg (-89:2:89), LonDeg (0:2:358)  comparison grid
%     Quantity ("geoid")   synthesis quantity; "ewh" etc. need kn
%     kn ([])              load Love numbers (user-supplied, 0..nmax)
%     Mask ([])            nlat x nlon logical region restriction
%     Names (["solution 1","solution 2"])
%     Plot (false)
%
%   Outputs
%     rep        (1,1) struct  fields:
%                  .spectral (1,1) struct  shx.diffSpectrum output
%                  .spatial  (1,1) struct  shx.spatialStats output
%                  .chi2dof  (1,1) double  error-realism statistic (NaN
%                                          when either solution lacks sigmas)
%                  .nmax     (1,1) double  common maximum degree used
%                  .quantity (1,1) string  compared quantity
%                  .names    (1,2) string  solution labels
%                  .rescaled (1,1) logical true when GM/R were unified
%     h          (1,1) graphics handle  figure handle (Plot = true only)
%
%   Example
%     rep = shx.compareSolutions(g, g.gaussian(350), Names = ["raw", "G350"]);
%     fprintf("crossover n = %d, spatial corr = %.3f\n", ...
%         rep.spectral.ncross, rep.spatial.corr)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.6.0).
arguments
    g1 (1,1) shCoefficients
    g2 (1,1) shCoefficients
    opts.LatDeg (1,:) double = -89:2:89
    opts.LonDeg (1,:) double = 0:2:358
    opts.Quantity (1,1) string = "geoid"
    opts.kn double = []
    opts.Mask = []
    opts.Names (1,2) string = ["solution 1", "solution 2"]
    opts.Plot (1,1) logical = false
end
% ---- common basis
nmax = min(g1.nmax, g2.nmax);
if g1.nmax > nmax, g1 = g1.truncate(nmax); end
if g2.nmax > nmax, g2 = g2.truncate(nmax); end
rescaled = false;
if abs(g1.GM - g2.GM) > 0 || abs(g1.R - g2.R) > 0
    g2 = g2.toReference(GM = g1.GM, R = g1.R);
    rescaled = true;
end
% ---- spectral
spec = shx.diffSpectrum(g1.C, g1.S, g2.C, g2.S);
% ---- error realism (coefficient domain)
chi2dof = NaN;
if ~isempty(g1.sigmaC) && ~isempty(g2.sigmaC)
    idx = shx.shIndex(nmax);
    d = shx.vecFromCS(g1.C - g2.C, g1.S - g2.S, idx);
    s2 = shx.vecFromCS(g1.sigmaC, g1.sigmaS, idx).^2 + ...
         shx.vecFromCS(g2.sigmaC, g2.sigmaS, idx).^2;
    use = isfinite(d) & isfinite(s2) & s2 > 0;
    if any(use), chi2dof = mean(d(use).^2 ./ s2(use)); end
end
% ---- spatial
A = g1.synthesis(opts.LatDeg, opts.LonDeg, quantity = opts.Quantity, ...
    kn = opts.kn);
B = g2.synthesis(opts.LatDeg, opts.LonDeg, quantity = opts.Quantity, ...
    kn = opts.kn);
st = shx.spatialStats(A, B, opts.LatDeg, opts.LonDeg, Mask = opts.Mask);
rep = struct('spectral', spec, 'spatial', st, 'chi2dof', chi2dof, ...
    'nmax', nmax, 'quantity', opts.Quantity, 'names', opts.Names, ...
    'rescaled', rescaled);
% ---- visualization
h = gobjects(0);
if opts.Plot
    h = figure('Color', 'w');
    tiledlayout(h, 2, 2, 'Padding', 'compact');
    nexttile;
    semilogy(spec.n, spec.amp1, 'k-', spec.n, spec.amp2, 'b-', ...
        spec.n, spec.diffAmp, 'r-', 'LineWidth', 1.1);
    grid on; xlabel('degree n'); ylabel('degree amplitude');
    legend([opts.Names, "difference"], 'Location', 'best');
    if isfinite(spec.ncross)
        xline(spec.ncross, ':', sprintf('n_c = %d', spec.ncross));
    end
    title('spectra');
    nexttile;
    plot(spec.n, spec.degCorr, 'k.-'); grid on; ylim([-0.1 1.05]);
    xlabel('degree n'); ylabel('degree correlation'); title('pattern per degree');
    nexttile;
    imagesc(opts.LonDeg, opts.LatDeg, A - B); axis xy; colorbar;
    xlabel('longitude'); ylabel('latitude');
    title(sprintf('%s difference (rmsd %.3g)', opts.Quantity, st.rmsd));
    nexttile;
    so = shx.shOrderRMS(g1.C - g2.C, g1.S - g2.S);
    semilogy(so.order, so.ordRMS, 'r.-'); grid on;
    xlabel('order m'); ylabel('order RMS of difference');
    title('difference per order (stripes/aliases)');
end
end
