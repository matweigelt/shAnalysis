function y = icgemDate2Year(d)
%ICGEMDATE2YEAR Convert ICGEM/GRACE yyyymmdd(.hhmm) epochs to decimal years.
%
%   Y = shLowLevel.icgemDate2Year(D) converts date codes of the form yyyymmdd or
%   yyyymmdd.hhmm (as used in ICGEM gfct files, TN-13/TN-14) to decimal
%   years, using day-of-year / (days in year) with proper leap years.
%   Values that do not look like date codes (< 1.8e7, e.g. already-decimal
%   years in legacy synthetic files) are returned UNCHANGED - this rule
%   keeps pre-v2.1 files evaluating identically.
%
%   Inputs   d  double array, yyyymmdd(.hhmm) or decimal years
%   Inputs
%     d               dates as yyyymmdd(.hhmm) numbers or decimal
%                     years (passed through)
%
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     y          (same size as input) double   decimal years
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    d double
end
y = d;
isDate = d >= 1.8e7;
if ~any(isDate(:)), return; end
dd = d(isDate);
di = floor(dd);
yyyy = floor(di / 1e4);
mm   = floor(mod(di, 1e4) / 1e2);
day  = mod(di, 1e2);
frac = dd - di;                                % .hhmm
hh   = floor(frac * 1e2);
mins = round(mod(frac * 1e4, 1e2));
% pure arithmetic (v3.1.1): the previous datetime/calendarDuration
% construction cost ~1 ms per call - fatal inside shReadGFC's line
% parser, where variable-term files (GRGS mean fields) call this for
% every record (tens of minutes for a 73 MB model). Same values.
cumdays = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
leap = (mod(yyyy, 4) == 0 & mod(yyyy, 100) ~= 0) | mod(yyyy, 400) == 0;
doy = reshape(cumdays(mm), size(mm)) + day - 1 + (mm > 2) .* leap ...
    + (hh + mins / 60) / 24;
y(isDate) = yyyy + doy ./ (365 + leap);
end
