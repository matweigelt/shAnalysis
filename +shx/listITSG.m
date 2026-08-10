function [T, info] = listITSG(opts)
%LISTITSG Catalogue of ITSG gravity field series on the TU Graz server.
%
%   T = shx.listITSG() enumerates the available ITSG releases and
%   product folders (releases ITSG-Grace2014/2016/2018 and the
%   continuously extended ITSG-Grace_operational; products monthly at
%   n60/n96/n120 and daily Kalman solutions) by parsing the server
%   directory indices - so mixing files from different releases, the
%   trap of guessing folder names, becomes impossible. The catalogue
%   row number is the selection handle for shx.fetchITSG(Catalog=...).
%   Behaviour mirrors shx.listICGEM.
%
%   Options
%     BaseURL ("https://ftp.tugraz.at/pub/ITSG/GRACE")  server base, or
%              a LOCAL MIRROR FOLDER (offline/testing: directory tree
%              is walked instead of scraped)
%     Timeout (30)
%
%   Outputs
%     T          (K x 5) table  columns:
%                  idx     (1,1) double  selection number for fetchITSG
%                  release (1,1) string  e.g. "ITSG-Grace2018"
%                  product (1,1) string  "monthly" | "daily" | "static"
%                  nmax    (1,1) double  60/96/120 (NaN for daily/static)
%                  url     (1,1) string  folder URL or mirror path
%     info       (1,1) struct  fields: source (string), nReleases (double)
%
%   Example
%     T = shx.listITSG();
%     disp(T)                            % pick rows, then:
%     % shx.fetchITSG(Catalog = [3 4])   % fetch those folders completely
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-10 (v3.0.0).
arguments
    opts.BaseURL (1,1) string = "https://ftp.tugraz.at/pub/ITSG/GRACE"
    opts.Timeout (1,1) double = 30
end
base = opts.BaseURL;
local = isfolder(base);
rows = {};
if local
    d = dir(char(base));
    rel = string({d([d.isdir]).name});
    rel = rel(startsWith(rel, "ITSG-Grace"));
else
    rel = hrefDirs(base + "/", opts.Timeout);
    rel = rel(startsWith(rel, "ITSG-Grace"));
end
for r = rel
    sub = subDirs(base, r, local, opts.Timeout);
    if any(sub == "monthly")
        mv = subDirs(base, r + "/monthly", local, opts.Timeout);
        for v = mv(startsWith(mv, "monthly_n"))
            nm = double(string(extractAfter(v, "monthly_n")));
            if isfinite(nm)
                rows(end+1, :) = {r, "monthly", nm, ...
                    join([base, r, "monthly", v], "/")}; %#ok<AGROW>
            end
        end
    end
    if any(sub == "daily_kalman")
        rows(end+1, :) = {r, "daily", NaN, ...
            join([base, r, "daily_kalman"], "/")}; %#ok<AGROW>
    end
    if any(sub == "static")
        rows(end+1, :) = {r, "static", NaN, ...
            join([base, r, "static"], "/")}; %#ok<AGROW>
    end
end
if isempty(rows)
    error('shx:listITSG:empty', 'No ITSG releases found under %s.', base);
end
K = size(rows, 1);
T = table((1:K)', string(rows(:, 1)), string(rows(:, 2)), ...
    cell2mat(rows(:, 3)), string(rows(:, 4)), ...
    'VariableNames', {'idx', 'release', 'product', 'nmax', 'url'});
info = struct('source', base, 'nReleases', numel(rel));
end

function out = subDirs(base, rel, local, timeout)
if local
    d = dir(fullfile(char(base), char(strrep(rel, "/", filesep))));
    out = string({d([d.isdir]).name});
    out = out(~startsWith(out, "."));
else
    out = hrefDirs(base + "/" + rel + "/", timeout);
end
end

function names = hrefDirs(url, timeout)
html = webread(url, weboptions('Timeout', timeout));
tok = regexp(html, 'href="([^"/?][^"?]*)/"', 'tokens');
names = string(cellfun(@(t) t{1}, tok, 'UniformOutput', false));
names = unique(names, 'stable');
end
