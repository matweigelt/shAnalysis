function out = makeTutorials(opts)
%MAKETUTORIALS Generate Live Script tutorials from the demo registry.
%
%   FILES = shLowLevel.makeTutorials() writes one tutorial per demo case
%   into tutorials/ as a plain .m script with Live Editor section
%   markers, and converts each to .mlx when the converter is available.
%
%   The tutorials are GENERATED, not hand-written. demo_shAnalysis is
%   the single source of truth for what the toolbox demonstrates; a
%   parallel set of hand-maintained .mlx files would drift from it
%   within one release, and .mlx is a binary zip, so the drift would not
%   even be visible in a diff. Regenerating is cheap; keeping two copies
%   in step by hand is not.
%
%   Each tutorial contains the demo's title, the functions it exercises,
%   the explanatory text from the case registry, and a runnable call.
%   They are meant as an entry point for someone who has just installed
%   the toolbox - "open this, press Run, then change something" - not as
%   a replacement for the workflow guide, which explains the theory.
%
%   Options
%     Dest (fullfile(shLowLevel.version().Root, "tutorials"))  output
%             folder; created if absent
%     Cases ("all")  which demo cases to write ("all", "core", or an
%             array of ids such as ["D01" "D05"])
%     Convert (true)  also write .mlx. The converter is an internal
%             MATLAB API and is not present in every release, so a
%             failure is reported and the .m files are kept: those are
%             the useful artefact either way, and they are the ones
%             worth putting under version control
%     Quiet (false)  suppress progress output
%
%   Outputs
%     out        (1,1) struct  fields: files (1,N string, the .m files
%                written), mlx (1,M string, the .mlx files written -
%                empty when Convert is false or unavailable), dest
%                (1,1 string), converted (1,1 logical)
%
%   Example
%     shLowLevel.makeTutorials(Cases = ["D01" "D05"]);
%     edit(fullfile("tutorials", "T_D01_read_and_spectra.m"))
%
%   See also demo_shAnalysis, shLowLevel.version.
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.8.0).
arguments
    opts.Dest (1,1) string = ""
    opts.Cases (1,:) string = "all"
    opts.Convert (1,1) logical = true
    opts.Quiet (1,1) logical = false
end
v = shLowLevel.version();
dest = opts.Dest;
if strlength(dest) == 0
    dest = fullfile(v.Root, "tutorials");
end
if ~isfolder(dest)
    mkdir(dest);
end
% "list" PRINTS the table as well as returning it, which would spam a
% generator run and every test; capture the output instead
evalc('reg = demo_shAnalysis("list");');
ids = string({reg.id});
sel = opts.Cases;
if isscalar(sel) && sel == "all"
    sel = ids;
elseif isscalar(sel) && sel == "core"
    sel = ["D01" "D02" "D04" "D05" "D15"];
end
unknown = setdiff(sel, ids);
assert(isempty(unknown), 'shLowLevel:makeTutorials:unknownCase', ...
    'Unknown demo case(s): %s. Use demo_shAnalysis("list").', ...
    strjoin(unknown, ", "));

files = strings(1, 0);
mlx = strings(1, 0);
converted = true;
for k = 1:numel(sel)
    r = reg(ids == sel(k));
    name = "T_" + r.id + "_" + slug(r.title);
    fm = fullfile(dest, name + ".m");
    writeTutorial(fm, r, v);
    files(end+1) = string(fm); %#ok<AGROW>
    if opts.Convert && converted
        fx = fullfile(dest, name + ".mlx");
        try
            matlab.internal.livecode.FileModel.convertFileToLiveCode( ...
                char(fm), char(fx));
            mlx(end+1) = string(fx); %#ok<AGROW>
        catch err
            % The converter is an internal API. If it is missing the .m
            % files are still the useful artefact, so say so once and
            % carry on rather than failing the whole run.
            converted = false;
            warning('shLowLevel:makeTutorials:noConverter', ...
                ['Live Script conversion unavailable (%s). The .m ' ...
                 'tutorials were written and can be opened in the Live ' ...
                 'Editor with "Save As" to produce .mlx by hand.'], ...
                err.identifier);
        end
    end
    if ~opts.Quiet
        fprintf('  %s  %s\n', r.id, name);
    end
end
if ~opts.Quiet
    fprintf('%d tutorial(s) in %s%s\n', numel(files), dest, ...
        stringIf(~isempty(mlx), sprintf(' (+%d .mlx)', numel(mlx))));
end
out = struct('files', files, 'mlx', mlx, 'dest', string(dest), ...
    'converted', converted && opts.Convert);
end

% ------------------------------------------------------------- helpers
function writeTutorial(path, r, v)
%WRITETUTORIAL One tutorial: headings, prose, a runnable call.
L = strings(0, 1);
L(end+1) = "%% " + r.id + " - " + r.title;
L(end+1) = "% A guided tour of one shAnalysis capability, generated from the";
L(end+1) = "% demo registry (demo_shAnalysis). Press Run, then change something.";
L(end+1) = "%";
L(end+1) = "% Functions exercised: " + r.fns;
L(end+1) = "";
L(end+1) = "%% Set up";
L(end+1) = "% setup_shAnalysis puts the toolbox on the path and, with";
L(end+1) = "% Download=, fetches the data files the demos use. Everything";
L(end+1) = "% below works offline against the shipped fixtures.";
L(end+1) = "setup_shAnalysis;";
L(end+1) = "";
L(end+1) = "%% Run the case";
L(end+1) = "% Visible=true keeps the figures on screen; the demo prints what";
L(end+1) = "% it is doing as it goes.";
L(end+1) = "demo_shAnalysis(""" + r.id + """, Visible = true);";
L(end+1) = "";
L(end+1) = "%% Where to go next";
L(end+1) = "% * |doc shAnalysis| for the topic pages on these functions";
L(end+1) = "% * The workflow guide (docs/) explains the theory and the";
L(end+1) = "%   conventions behind them, with worked recipes.";
L(end+1) = "% * |demo_shAnalysis(""list"")| lists every case.";
L(end+1) = "%";
L(end+1) = "% Generated by shLowLevel.makeTutorials for shAnalysis " + ...
    v.Version + ". Edits here are";
L(end+1) = "% overwritten on the next run - change demo_shAnalysis instead.";
L(end+1) = "%";
L(end+1) = "% Developed by Matthias Weigelt with the help of Claude.";
fid = fopen(path, 'w');
assert(fid > 0, 'shLowLevel:makeTutorials:cannotWrite', ...
    'Cannot write %s', path);
cl = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', L{:});
end

function s = slug(t)
%SLUG A safe, short file-name fragment from a title.
%   repmat("x", 1, 40) builds a 1x40 STRING ARRAY, not a 40-character
%   string - the padding trick that works for char silently produces an
%   array here, and the result then cannot be a file name. Truncate
%   directly instead.
s = lower(string(t));
s = regexprep(s, '[^a-z0-9]+', '_');
s = regexprep(s, '^_+|_+$', '');
if strlength(s) > 40
    s = extractBefore(s, 41);
end
if strlength(s) == 0
    s = "case";
end
end

function s = stringIf(tf, t)
if tf, s = t; else, s = ""; end
end
