function [ts, rep] = standardChain(folder, opts)
%STANDARDCHAIN Canonical GRACE/GRACE-FO post-processing in one call.
%
%   [TS, REP] = shLowLevel.standardChain(FOLDER) encodes the standard
%   chain as a single, correctly ordered entry point:
%     1. read all monthly solutions from FOLDER (shSeries.fromFolder)
%     2. replace C20/C30 with the SLR values (TN-14)
%     3. add the degree-1 (geocenter) coefficients (TN-13)
%     4. optionally subtract a GIA trend field
%     5. filter (Gaussian radius, a DDK, or a shLowLevel.designFilter W)
%   The order matters (TN corrections BEFORE filtering, GIA on the
%   corrected series) and getting it wrong is the classic beginner
%   error this function exists to prevent. Every step is recorded in
%   REP and in the series history; steps are skippable for
%   sensitivity studies. Data files come from the persistent data
%   folder by default (setup_shAnalysis fetches them) and can be
%   overridden with explicit paths for reproducible pipelines.
%
%   Inputs
%     folder  (1,1) string  folder of monthly .gfc(.gz) files
%
%   Options
%     Pattern ("*.gfc*")   file pattern for shSeries.fromFolder
%     Truncate (NaN)       truncate to this nmax after reading
%     TN14 (true)          apply SLR C20/C30
%     TN14File ("")        explicit TN-14 path (default: dataFolder copy)
%     Degree1 ("GFZ")      TN-13 provider "GFZ" | "CSR" | "JPL" | "none"
%     Degree1File ("")     explicit TN-13 path (default: dataFolder copy)
%     GIA ([])             shLowLevel trend field (shCoefficients, [1/yr])
%                          subtracted as (t - GIAEpoch) * GIA
%     GIAEpoch (NaN)       reference epoch of the GIA trend (NaN: mean
%                          epoch of the series)
%     Filter ("gauss300")  "none" | "gaussN" (N = radius [km]) | "DDKn"
%                          (n = 1..8, fetched on demand) | a W struct
%                          from shLowLevel.readDDK / designFilter
%     Tolerance (0.05)     maximum |epoch difference| [yr] accepted when
%                          matching a TN-13/TN-14 record to a solution
%     OnMissing ("drop")   what to do with solution epochs the correction
%                          tables do not reach. The tables ALWAYS trail
%                          the solutions - a provider publishes TN-13
%                          weeks after the monthly field - so the newest
%                          months of a fresh series are routinely
%                          uncovered and the old behaviour was to stop
%                          with an error. "drop" removes them and records
%                          it in REP; "error" restores the old behaviour.
%                          Dropping is the default because the only other
%                          option that keeps them would splice
%                          UNCORRECTED months onto corrected ones, which
%                          puts a step in the series exactly where people
%                          look for the newest signal
%     Quiet (false)        suppress progress output
%
%   Outputs
%     ts         (1,1) shSeries  the processed series (full history)
%     rep        (1,1) struct  fields: files (1,T string), steps
%                (1,K string, applied stages with parameters), nmax
%                (1,1 double), epochs (T,1 double), version (1,1 string,
%                shLowLevel.version at run time), created (1,1 string)
%
%   Example
%     [ts, rep] = shLowLevel.standardChain("D:/grace/monthly", ...
%         Filter = "DDK3", Degree1 = "GFZ");
%     disp(rep.steps')
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-10 (v3.1.0).
arguments
    folder (1,1) string
    opts.Pattern (1,1) string = "*.gfc*"
    opts.Truncate (1,1) double = NaN
    opts.TN14 (1,1) logical = true
    opts.TN14File (1,1) string = ""
    opts.Degree1 (1,1) string {mustBeMember(opts.Degree1, ...
        ["GFZ", "CSR", "JPL", "none"])} = "GFZ"
    opts.Degree1File (1,1) string = ""
    opts.GIA = []
    opts.GIAEpoch (1,1) double = NaN
    opts.Filter = "gauss300"
    opts.Tolerance (1,1) double {mustBePositive} = 0.05
    opts.OnMissing (1,1) string {mustBeMember(opts.OnMissing, ...
        ["drop", "error"])} = "drop"
    opts.Quiet (1,1) logical = false
end
steps = strings(1, 0);
say = @(s) fprintf('  %s\n', s);
if opts.Quiet, say = @(s) []; end %#ok<NASGU>
% ---- 1. read
ts = shSeries.fromFolder(folder, Pattern = opts.Pattern);
if isfinite(opts.Truncate), ts = ts.truncate(opts.Truncate); end
files = strings(1, 0);
d = dir(fullfile(char(folder), char(opts.Pattern)));
if ~isempty(d), files = string({d.name}); end
steps(end+1) = sprintf("read %d epochs (nmax %d)", ts.nEpochs, ts.nmax);
if ~opts.Quiet, fprintf('  %s\n', steps(end)); end
% ---- 2. TN-14
if opts.TN14
    tnf = opts.TN14File;
    if strlength(tnf) == 0
        tnf = fullfile(shLowLevel.dataFolder(), 'tn', ...
            'TN-14_C30_C20_SLR_GSFC.txt');
    end
    if ~isfile(tnf)
        error('shLowLevel:standardChain:noTN14', ...
            ['TN-14 file not found (%s). Run setup_shAnalysis or give ' ...
             'TN14File= explicitly.'], tnf);
    end
    tn = shLowLevel.readTN14(tnf);
    [ts, nDrop] = dropUncovered(ts, tn.epoch, opts.Tolerance, ...
        opts.OnMissing, "TN-14");
    if nDrop > 0
        steps(end+1) = sprintf("dropped %d epoch(s) beyond TN-14 coverage (table ends %.3f)", nDrop, max(tn.epoch));
        if ~opts.Quiet, fprintf('  %s\n', steps(end)); end
    end
    ts = ts.applyTN14(tn, Tolerance = opts.Tolerance);
    steps(end+1) = "TN-14 C20/C30 (" + string(tnf) + ")";
    if ~opts.Quiet, fprintf('  %s\n', steps(end)); end
end
% ---- 3. degree 1
if opts.Degree1 ~= "none"
    d1f = opts.Degree1File;
    if strlength(d1f) == 0
        d1f = fullfile(shLowLevel.dataFolder(), 'tn', ...
            sprintf('TN-13_GEOC_%s_RL06.3.txt', opts.Degree1));
    end
    if ~isfile(d1f)
        error('shLowLevel:standardChain:noTN13', ...
            ['TN-13 file not found (%s). Run setup_shAnalysis or give ' ...
             'Degree1File= explicitly.'], d1f);
    end
    tn13 = shLowLevel.readTN13(d1f);
    [ts, nDrop] = dropUncovered(ts, tn13.epoch, opts.Tolerance, ...
        opts.OnMissing, "TN-13");
    if nDrop > 0
        steps(end+1) = sprintf("dropped %d epoch(s) beyond TN-13 coverage (table ends %.3f)", nDrop, max(tn13.epoch));
        if ~opts.Quiet, fprintf('  %s\n', steps(end)); end
    end
    ts = ts.addDegree1(tn13, Tolerance = opts.Tolerance);
    steps(end+1) = "degree-1 " + opts.Degree1 + " (" + string(d1f) + ")";
    if ~opts.Quiet, fprintf('  %s\n', steps(end)); end
end
% ---- 4. GIA
if ~isempty(opts.GIA)
    if ~isa(opts.GIA, 'shCoefficients')
        error('shLowLevel:standardChain:badGIA', ...
            'GIA must be an shCoefficients trend field [1/yr].');
    end
    t0 = opts.GIAEpoch;
    if ~isfinite(t0), t0 = mean(ts.epochs); end
    gT = opts.GIA;
    if gT.nmax > ts.nmax, gT = gT.truncate(ts.nmax); end
    n1 = ts.nmax + 1; T = ts.nEpochs;
    Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
    for k = 1:T
        dt = ts.epochs(k) - t0;
        Cs(:, :, k) = dt * gT.C; Ss(:, :, k) = dt * gT.S;
    end
    ts = ts - shSeries(Cs, Ss = Ss, Epochs = ts.epochs);
    steps(end+1) = sprintf("GIA trend removed (t0 = %.2f)", t0);
    if ~opts.Quiet, fprintf('  %s\n', steps(end)); end
end
% ---- 5. filter
F = opts.Filter;
if isstruct(F)
    ts = ts.applyDDK(F);
    steps(end+1) = "filter: custom W (" + string(F.name) + ")";
elseif isstring(F) || ischar(F)
    F = string(F);
    if F == "none"
        steps(end+1) = "filter: none";
    elseif startsWith(F, "gauss")
        r = double(extractAfter(F, "gauss"));
        assert(isfinite(r) && r > 0, 'shLowLevel:standardChain:badFilter', ...
            'Gaussian filter must be "gauss<radiusKm>" (got "%s").', F);
        ts = ts.gaussian(r);
        steps(end+1) = sprintf("filter: Gaussian %g km", r);
    elseif startsWith(F, "DDK")
        Wf = shLowLevel.readDDK(char(F));
        ts = ts.applyDDK(Wf);
        steps(end+1) = "filter: " + F;
    else
        error('shLowLevel:standardChain:badFilter', ...
            'Filter must be "none", "gauss<km>", "DDKn", or a W struct.');
    end
else
    error('shLowLevel:standardChain:badFilter', ...
        'Filter must be a string or a W struct.');
end
if ~opts.Quiet, fprintf('  %s\n', steps(end)); end
v = shLowLevel.version();
rep = struct('files', files, 'steps', steps, 'nmax', ts.nmax, ...
    'epochs', ts.epochs, 'version', string(v.Name + " " + v.Version), ...
    'created', string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
end

% ------------------------------------------------------------- helpers
function [ts, nDropped] = dropUncovered(ts, tableEpochs, tol, mode, what)
%DROPUNCOVERED Remove solution epochs no correction record can reach.
lo = min(tableEpochs) - tol;
hi = max(tableEpochs) + tol;
keep = ts.epochs >= lo & ts.epochs <= hi;
nDropped = nnz(~keep);
if nDropped == 0
    return
end
if mode == "error"
    error('shLowLevel:standardChain:uncoveredEpochs', ...
        ['%d solution epoch(s) lie outside the %s table (%.3f..%.3f). ' ...
         'The correction tables trail the solutions, so this is normal ' ...
         'for a fresh series: OnMissing="drop" (the default) excludes ' ...
         'them, or supply an updated table.'], ...
        nDropped, what, min(tableEpochs), max(tableEpochs));
end
if ~any(keep)
    error('shLowLevel:standardChain:noCoveredEpochs', ...
        ['No solution epoch falls inside the %s table (%.3f..%.3f) - ' ...
         'the table and the series do not overlap at all.'], ...
        what, min(tableEpochs), max(tableEpochs));
end
ts = ts.select(keep);
end
