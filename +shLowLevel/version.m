function v = version()
%VERSION Toolbox metadata: name, version, date, root path, provenance.
%
%   V = shLowLevel.version() returns the toolbox metadata as a struct. The
%   version string and date are PARSED FROM Contents.m (the MATLAB
%   toolbox convention line "% Version x.y.z ... dd-Mmm-yyyy"), so
%   there is a single source of truth also honoured by MATLAB's own
%   ver('shAnalysis'). Called without an output, prints a summary.
%
%   Outputs
%     v          (1,1) struct  fields:
%                  .Name       (1,1) string  "shAnalysis"
%                  .Version    (1,1) string  e.g. "2.5.1"
%                  .Date       (1,1) string  release date from Contents.m
%                  .Root       (1,1) string  toolbox root folder
%                  .Provenance (1,1) string  authorship line
%
%   Example
%     v = shLowLevel.version();
%     fprintf("%s %s (%s)\n", v.Name, v.Version, v.Date)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.5.1).
root = fileparts(fileparts(mfilename('fullpath')));
cf = fullfile(root, 'Contents.m');
verStr = "unknown"; datStr = "unknown";
if isfile(cf)
    txt = fileread(cf);
    tok = regexp(txt, '^%\s*Version\s+(\S+).*?(\d{2}-\w{3}-\d{4})\s*$', ...
        'tokens', 'once', 'lineanchors');
    if ~isempty(tok)
        verStr = string(tok{1}); datStr = string(tok{2});
    end
end
out = struct('Name', "shAnalysis", 'Version', verStr, 'Date', datStr, ...
    'Root', string(root), 'Provenance', ...
    "Developed by Matthias Weigelt with the help of Claude (Fable 5)");
if nargout == 0
    fprintf('%s %s (%s)\n  root: %s\n  %s\n', out.Name, out.Version, ...
        out.Date, out.Root, out.Provenance);
else
    v = out;
end
end
