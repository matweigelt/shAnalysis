function tn = readTN14(filename)
%READTN14 Parse a TN-14 style SLR C20/C30 replacement table.
%
%   TN = shLowLevel.readTN14(FILENAME) reads the CSR/GFZ "TN-14" technical-note
%   text format: an arbitrary text header followed by data lines with 10
%   whitespace-separated numeric columns
%       MJD_begin  yearfrac_begin  C20  C20-mean*1e10  sigC20*1e10 ...
%       C30  C30-mean*1e10  sigC30*1e10  MJD_end  yearfrac_end
%   C30 columns contain 'NaN' for months where no SLR C30 is provided
%   (pre GRACE-FO accelerometer degradation).
%
%   ASSUMPTION: the 10-column layout above. The parser skips any line
%   that does not yield exactly 10 numeric tokens (NaN accepted), so
%   headers of any length are tolerated; verify tn.epoch range against
%   the file's own documentation after reading unfamiliar versions.
%
%   Inputs
%     filename  char/string
%   Outputs (struct TN), K = number of data lines
%     epoch      (K,1) double   mid-epoch, decimal years
%     C20        (K,1) double
%     sigmaC20   (K,1) double
%     C30        (K,1) double   NaN where not provided
%     sigmaC30   (K,1) double   NaN where not provided
%     epochStart, epochStop (K,1) double
%
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     tn         struct: epoch/epochStart/epochStop (M x 1 decimal years), C20, sigmaC20, C30, sigmaC30 (M x 1; C30 NaN before 2012.16)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    filename {mustBeTextScalar}
end

fid = fopen(filename, 'r');
assert(fid > 0, 'shLowLevel:readTN14:fileNotFound', 'Cannot open %s.', char(filename));
cl = onCleanup(@() fclose(fid));

rows = [];
line = fgetl(fid);
while ischar(line)
    v = sscanf(strrep(line, 'NaN', 'nan'), '%f');
    if numel(v) == 10
        rows(end+1, :) = v'; %#ok<AGROW>
    end
    line = fgetl(fid);
end
assert(~isempty(rows), 'shLowLevel:readTN14:noData', ...
    'No 10-column numeric data lines found in %s.', char(filename));

tn.epochStart = rows(:, 2);
tn.epochStop  = rows(:, 10);
tn.epoch      = (tn.epochStart + tn.epochStop) / 2;
tn.C20        = rows(:, 3);
tn.sigmaC20   = rows(:, 5)  * 1e-10;
tn.C30        = rows(:, 6);
tn.sigmaC30   = rows(:, 8)  * 1e-10;
end
