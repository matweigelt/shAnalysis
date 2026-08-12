function p = gravisDataFolder()
%GRAVISDATAFOLDER Path of the shipped GravIS auxiliary data.
%
%   P = shLowLevel.gravisDataFolder() returns the absolute path of the
%   data/gravis folder shipped with the toolbox: frozen (2026-08-12)
%   copies of the GravIS COST-G Level-2B auxiliary files and the GravIS
%   basin polygons that gravisL2B, greenlandChain, antarcticaChain and
%   twsChain use when no gravisFolder is given. The correction tables
%   trail the monthly solutions - for epochs after the freeze fetch
%   fresh copies (see data/gravis/README.md) and pass gravisFolder=.
%
%   Inputs
%     (none)
%   Outputs
%     p  (1,1) string  absolute folder path
%
%   Example
%     dir(shLowLevel.gravisDataFolder())
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.8.8).
here = fileparts(mfilename('fullpath'));           % .../+shLowLevel
p = string(fullfile(fileparts(here), 'data', 'gravis'));
if ~isfolder(p)
    error('shLowLevel:gravisDataFolder:missing', ...
        'Shipped data folder not found: %s (broken installation?).', p);
end
end
