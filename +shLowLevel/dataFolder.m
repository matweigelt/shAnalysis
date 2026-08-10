function folder = dataFolder(newFolder)
%DATAFOLDER Get or set the persistent user data folder of the toolbox.
%
%   F = shLowLevel.dataFolder()          returns the current data folder
%                                 (created on demand). Default:
%                                 <toolbox root>/data.
%   F = shLowLevel.dataFolder(PATH)      sets and persists a user folder
%                                 (MATLAB preference, survives sessions:
%                                 setpref('shAnalysis','dataFolder')).
%   F = shLowLevel.dataFolder("reset")   back to the default.
%
%   Consumers: shLowLevel.fetchITSG (itsg_series/), shLowLevel.fetchDDK (DDK/),
%   shLowLevel.fetchICGEM (icgem/), shLowLevel.readDDK("DDK<n>") name resolution, and
%   the demo's real-data discovery. Point this at a shared network
%   location to reuse downloads across machines.
%
%   Claude (Fable 5), 2026-08-07 (v2.4.1).
%
%   Inputs
%     newFolder  char/string  new persistent data folder (omit to query the current one)
%   Outputs
%     folder     (1 x 1) string   current data folder (created on demand)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    newFolder (1,1) string = ""
end
group = 'shAnalysis'; prefName = 'dataFolder';
if newFolder == "reset"
    if ispref(group, prefName), rmpref(group, prefName); end
    newFolder = "";
end
if strlength(newFolder) > 0
    setpref(group, prefName, char(newFolder));
end
if ispref(group, prefName)
    folder = string(getpref(group, prefName));
else
    root = fileparts(fileparts(mfilename('fullpath')));   % .../+shLowLevel -> root
    folder = string(fullfile(root, 'data'));
end
if ~isfolder(folder), mkdir(folder); end
end
