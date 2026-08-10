function writeGFC(filename, C, S, GM, R, opts)
%WRITEGFC Write Stokes coefficients to an ICGEM ASCII .gfc file.
%
%   shx.writeGFC(FILENAME, C, S, GM, R) writes a static ICGEM-format
%   gravity field file readable by shx.shReadGFC (round-trip tested to
%   full double precision: coefficients are written with %22.14e).
%
%   Inputs
%     filename   char/string
%     C, S       (nmax+1)x(nmax+1) double   fully normalized coefficients
%     GM         (1,1) double  [m^3/s^2]
%     R          (1,1) double  [m]
%     opts.SigmaC, opts.SigmaS  same layout, formal errors ([] -> zeros,
%                               errors flag written as 'no')
%     opts.ModelName   string   default "shAnalysis_export"
%     opts.TideSystem  string   default "unknown"
%     opts.Comment     string   extra header comment lines ("" = none)
%   Outputs: none (file on disk)
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    filename {mustBeTextScalar}
    C double
    S double
    GM (1,1) double {mustBePositive}
    R (1,1) double {mustBePositive}
    opts.SigmaC double = []
    opts.SigmaS double = []
    opts.ModelName {mustBeTextScalar} = "shAnalysis_export"
    opts.TideSystem {mustBeTextScalar} = "unknown"
    opts.Comment {mustBeTextScalar} = ""
end

nmax = size(C, 1) - 1;
assert(isequal(size(C), size(S)) && size(C,1) == size(C,2), ...
    'shx:writeGFC:badInput', 'C and S must be square and equally sized.');

hasSigma = ~isempty(opts.SigmaC) && ~isempty(opts.SigmaS);
sC = opts.SigmaC; sS = opts.SigmaS;
if ~hasSigma, sC = zeros(size(C)); sS = zeros(size(S)); end
sC(isnan(sC)) = 0; sS(isnan(sS)) = 0;

fid = fopen(filename, 'w');
assert(fid > 0, 'shx:writeGFC:cannotOpen', 'Cannot open %s for writing.', ...
    char(filename));
cl = onCleanup(@() fclose(fid));

if strlength(opts.Comment) > 0
    fprintf(fid, '# %s\n', char(opts.Comment));
end
fprintf(fid, '# written by shAnalysis (Claude, Fable 5), %s\n', ...
    char(datetime('now', Format = 'yyyy-MM-dd HH:mm')));
fprintf(fid, 'product_type            gravity_field\n');
fprintf(fid, 'modelname               %s\n', char(opts.ModelName));
fprintf(fid, 'earth_gravity_constant  %.10e\n', GM);
fprintf(fid, 'radius                  %.10e\n', R);
fprintf(fid, 'max_degree              %d\n', nmax);
fprintf(fid, 'errors                  %s\n', char(string(ternary(hasSigma, "formal", "no"))));
fprintf(fid, 'norm                    fully_normalized\n');
fprintf(fid, 'tide_system             %s\n', char(opts.TideSystem));
fprintf(fid, 'end_of_head =========================================\n');

for n = 0:nmax
    for m = 0:n
        fprintf(fid, 'gfc %5d %5d %22.14e %22.14e %12.4e %12.4e\n', ...
            n, m, C(n+1,m+1), S(n+1,m+1), sC(n+1,m+1), sS(n+1,m+1));
    end
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
