function varargout = legendreALF(varargin)
%LEGENDREALF v1 compatibility wrapper -> shx.legendreALF (unchanged behavior).
%   Kept so existing scripts and the legacy test suite run as-is; new code
%   should use the shCoefficients / shSeries classes.
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
[varargout{1:nargout}] = shx.legendreALF(varargin{:});
end
