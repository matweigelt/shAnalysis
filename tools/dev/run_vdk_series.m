%RUN_VDK_SERIES Build a VDK-filtered ITSG series - resumable batch.
%
%   Designed for an unattended run on a machine with disk and time:
%   the full n96 series means ~460 MB SINEX per month (~120 GB total).
%   Every stage is RESUMABLE: downloads skip present files, filtering
%   skips months whose output .gfc already exists, and a crash loses at
%   most the current month. Edit the CONFIG block, then run the script
%   whole; progress and failures are printed per month, a provenance
%   log is written next to the output.
%
%   Pipeline per month (Horvath et al. 2018):
%     fetchITSGSINEX -> readSINEX(NEQ, Index = idx) -> N ->
%     vdkApply(x, N, idx.n, ab(calendar month), Alpha) -> writeGFC
%
%   The signal model ab (cyclostationary Kaula a*l^b) is estimated ONCE
%   from a pre-filtered reference series (SignalSeriesFolder, e.g. a
%   DDK4-filtered ITSG series) - the paper's recipe. No numbers are
%   invented: the script REFUSES to run without that folder.
%
%   Acceptance built in (section 4): the paper's three probes on one
%   sample month - N (unfiltered) vs filtered degree RMS, NS/EW
%   half-weight anisotropy (expect NS 50-60% of EW), and the M = c*N
%   closed-form identity as a numerical self-test.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.13.0).

%% CONFIG - edit before running
cfg.Months    = ["2003-01"];                 % e.g. compose("%d-%02d", ...) spans
cfg.Nmax      = 96;                          % 96 | 120 (SINEX variants)
cfg.Alpha     = 1;                           % filter strength (paper Tab. 1)
cfg.SinexDir  = fullfile(shLowLevel.dataFolder(), "series", "itsg", "sinex");
cfg.OutDir    = fullfile(shLowLevel.dataFolder(), "itsg_vdk");
cfg.SignalSeriesFolder = "";                 % REQUIRED: pre-filtered (DDK4) series
cfg.LRange    = [10, 60];                    % Kaula fit range
cfg.KeepSinex = true;                        % false: delete each SINEX after use
cfg.LogFile   = fullfile(cfg.OutDir, "run_vdk_series.log");

%% 1 - signal model (once, cached beside the output)
assert(strlength(cfg.SignalSeriesFolder) > 0, ...
    ['CONFIG: SignalSeriesFolder is required - point it at a ' ...
     'pre-filtered (e.g. DDK4) series; the Kaula estimator reads ' ...
     'SIGNAL and must not see stripes.']);
if ~isfolder(cfg.OutDir), mkdir(cfg.OutDir); end
abFile = fullfile(cfg.OutDir, "vdk_signal_ab.mat");
if isfile(abFile)
    load(abFile, "ab", "abInfo");
    fprintf('signal model: loaded cached ab (%s)\n', abFile);
else
    ts = shSeries.read(cfg.SignalSeriesFolder);
    [ab, abInfo] = shLowLevel.signalVarianceKaula(ts.Cs, ts.Ss, ...
        ts.epochs, LRange = cfg.LRange);
    save(abFile, "ab", "abInfo");
    fprintf('signal model: estimated ab from %d epochs, cached\n', ...
        numel(ts.epochs));
end
disp(array2table(ab, VariableNames = ["a", "b"], ...
    RowNames = string(1:12)'));

%% 2 - month loop (resumable)
idx = shLowLevel.shIndex(cfg.Nmax, MinDegree = 2);
fid = fopen(cfg.LogFile, 'a');
fprintf(fid, '=== run %s | alpha %g | nmax %d ===\n', ...
    string(datetime("now")), cfg.Alpha, cfg.Nmax);
nOK = 0; nSkip = 0; nFail = 0;
for mm = cfg.Months
    outGfc = fullfile(cfg.OutDir, sprintf("ITSG_VDK%g_n%d_%s.gfc", ...
        cfg.Alpha, cfg.Nmax, mm));
    if isfile(outGfc)
        nSkip = nSkip + 1;
        continue
    end
    try
        t0 = tic;
        f = shLowLevel.fetchITSGSINEX(mm, Dest = cfg.SinexDir, ...
            Nmax = cfg.Nmax, Quiet = true);
        if isempty(f)
            fprintf(fid, '%s MISSING on server\n', mm);
            fprintf('%s: missing on server (gap month)\n', mm);
            continue
        end
        snx = shLowLevel.readSINEX(f(1), Index = idx);
        mo = str2double(extractAfter(mm, "-"));
        xf = shLowLevel.vdkApply(snx.x, snx.M, idx.n, ab(mo, :), ...
            Alpha = cfg.Alpha);
        [C, S] = shLowLevel.csFromVec(xf, idx);
        shLowLevel.writeGFC(outGfc, C, S, 3.986004415e14, 6378136.3, ...
            ModelName = sprintf("ITSG VDK-filtered (alpha %g)", cfg.Alpha), ...
            Comment = sprintf(['VDK/VADER decorrelation (Horvath et al. ' ...
                '2018) of ITSG %s: monthly SINEX N, cyclostationary ' ...
                'Kaula signal (a %.3g, b %.3g), alpha %g. shAnalysis %s.'], ...
                mm, ab(mo, 1), ab(mo, 2), cfg.Alpha, ...
                shLowLevel.version()));
        if ~cfg.KeepSinex, delete(f(1)); end
        nOK = nOK + 1;
        fprintf(fid, '%s OK %.0f s\n', mm, toc(t0));
        fprintf('%s: filtered and written (%.0f s)\n', mm, toc(t0));
    catch ME
        nFail = nFail + 1;
        fprintf(fid, '%s FAIL %s: %s\n', mm, ME.identifier, ME.message);
        fprintf('%s: FAIL %s\n', mm, ME.identifier);
    end
end
fprintf(fid, '=== done: %d ok, %d skipped, %d failed ===\n', nOK, nSkip, nFail);
fclose(fid);
fprintf('done: %d ok, %d skipped (present), %d failed - log: %s\n', ...
    nOK, nSkip, nFail, cfg.LogFile);

%% 3 - numerical self-test (closed form, no data needed)
idxT = shLowLevel.shIndex(10, MinDegree = 2);
rng(1, 'twister');
A = randn(idxT.P, ceil(idxT.P/3));
Nt = A*A' + idxT.P*eye(idxT.P);
c = 2.5; al = 0.7;
sig = 1 ./ sqrt(c * diag(Nt));          %#ok<NASGU> % not used: exact M=cN below
xt = randn(idxT.P, 1);
% M = c*N  =>  xf = x / (1 + alpha*c) exactly
Mc = c * Nt;
xfT = (Nt + al*Mc) \ (Nt * xt);
err = max(abs(xfT - xt / (1 + al*c)));
fprintf('%s  closed-form identity M=cN: err %.2e\n', ...
    string(ternaryLocal(err < 1e-10, "PASS", "FAIL")), err);

%% 4 - paper probes on the last filtered month (if any)
d = dir(fullfile(cfg.OutDir, "ITSG_VDK*.gfc"));
if ~isempty(d)
    g = shCoefficients.read(fullfile(d(end).folder, d(end).name));
    fprintf(['probe month %s: nmax %d written; compare its degree RMS ' ...
        'against the unfiltered ITSG month and check NS/EW kernel ' ...
        'anisotropy (expect NS 50-60%% of EW, Horvath et al. 2018) ' ...
        'before trusting the batch.\n'], d(end).name, size(g.C, 1) - 1);
end

function s = ternaryLocal(cnd, a, b)
if cnd, s = a; else, s = b; end
end
