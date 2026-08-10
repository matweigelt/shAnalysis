function [T, info] = listICGEM(opts)
%LISTICGEM List gravity field models available at ICGEM (GFZ).
%
%   T = shx.listICGEM() returns a table of the STATIC global models
%   from https://icgem.gfz.de/tom_longtime (name, year, max degree,
%   data sources, direct .gfc download URL) - feed a row's name into
%   shx.fetchICGEM to download.
%
%   [T, INFO] = shx.listICGEM(Type="temporal") returns the FULL temporal
%   catalogue parsed from the series-page links of /sl/temporal: every
%   GRACE/GRACE-FO/SLR/Swarm series of every center (~70+ series across
%   the groups 01_GRACE, 02_COST-G, 03_other, 04_SLR). Columns:
%     group    catalogue group (e.g. "01_GRACE", "03_other")
%     center   processing center (CSR, GFZ, JPL, COST-G, ITSG, ...)
%     series   series name incl. sub-listing (e.g. "ITSG-Grace2018/monthly")
%     path     the catalogue path "group/center/series" - feed this into
%              the Series option to list the individual files
%     url      the series page (human-readable listing)
%     zip      whole-series ZIP download (can be hundreds of MB!)
%
%   T = shx.listICGEM(Type="temporal", Series=PATH) fetches the series
%   page for PATH (a 'path' value from the catalogue table) and returns
%   the INDIVIDUAL files of that series (name, url) - single monthly
%   .gfc files are directly downloadable (websave), no whole-series ZIP
%   needed. Example:
%       T = shx.listICGEM(Type="temporal");
%       F = shx.listICGEM(Type="temporal", ...
%               Series="03_other/ITSG/ITSG-Grace2018/monthly");
%       websave("m1.gfc", F.url(1));
%   (For ITSG monthlies shx.fetchITSG remains the convenient route.)
%
%   Options
%     Type ("static")   "static" | "temporal"
%     Series ("")       temporal only: catalogue path of one series to
%                       list its individual files (see above)
%     Source ("")       local HTML file instead of the live page
%                       (offline/testing; the parsers are fixture-tested)
%     Timeout (60)
%   Outputs
%     T     table:  static:          name, year, degree, data, url
%                   temporal:        group, center, series, path, url, zip
%                   temporal+Series: name, url
%     info  struct: source, retrieved, note
%
%   v2.4.2: the temporal branch previously matched only "getseries" links
%   in the static HTML (release notes of 3 centers -> 3 rows) and missed
%   the /sp/ series-page catalogue entirely; single files ARE statically
%   listed on the series pages (the old "JS-only" caveat is obsolete).
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-07 (v2.4.1; temporal catalogue + Series listing v2.4.2).
arguments
    opts.Type (1,1) string {mustBeMember(opts.Type, ["static","temporal"])} = "static"
    opts.Series (1,1) string = ""
    opts.Source (1,1) string = ""
    opts.Timeout (1,1) double = 60
end
base = "https://icgem.gfz.de";
listFiles = strlength(opts.Series) > 0;
if listFiles && opts.Type ~= "temporal"
    error('shx:listICGEM:seriesNeedsTemporal', ...
        'The Series option applies to Type="temporal" only.');
end
if strlength(opts.Source) > 0
    html = fileread(opts.Source);
    src = opts.Source;
else
    page = base + "/tom_longtime";
    if opts.Type == "temporal", page = base + "/sl/temporal"; end
    if listFiles
        page = base + "/sp/" + strrep(opts.Series, " ", "%20");
    end
    html = webread(page, weboptions('Timeout', opts.Timeout, ...
        'ContentType', 'text', 'CharacterEncoding', 'UTF-8'));
    src = page;
end
info = struct('source', src, 'retrieved', datetime('now'), 'note', "");
if opts.Type == "static"
    rows = regexp(html, '<tr[^>]*>(.*?)</tr>', 'tokens');
    name = strings(0, 1); year = zeros(0, 1); degree = strings(0, 1);
    data = strings(0, 1); url = strings(0, 1);
    for k = 1:numel(rows)
        cells = regexp(rows{k}{1}, '<td[^>]*>(.*?)</td>', 'tokens');
        if numel(cells) < 8, continue; end
        u = regexp(rows{k}{1}, '(/getmodel/gfc/[^"]+\.gfc)', 'tokens', 'once');
        if isempty(u), continue; end
        clean = @(c) strtrim(regexprep(regexprep(c, '<[^>]+>', ' '), '\s+', ' '));
        name(end+1, 1) = string(clean(cells{3}{1})); %#ok<AGROW>
        year(end+1, 1) = str2double(clean(cells{4}{1})); %#ok<AGROW>
        degree(end+1, 1) = string(clean(cells{5}{1})); %#ok<AGROW>
        data(end+1, 1) = string(clean(cells{6}{1})); %#ok<AGROW>
        url(end+1, 1) = base + string(u{1}); %#ok<AGROW>
    end
    assert(~isempty(name), 'shx:listICGEM:parseFailed', ...
        'No models parsed from %s - page layout may have changed.', src);
    T = table(name, year, degree, data, url);
elseif listFiles
    % individual files of ONE series page: /getseries/<path>/<file>
    m = regexp(html, 'href="(/getseries/[^"]+\.(?:gfc|gz))"', 'tokens');
    name = strings(0, 1); url = strings(0, 1); seen = strings(0, 1);
    for k = 1:numel(m)
        r = string(m{k}{1});
        if any(seen == r), continue; end
        seen(end+1, 1) = r; %#ok<AGROW>
        parts = split(r, "/");
        name(end+1, 1) = parts(end); %#ok<AGROW>
        url(end+1, 1) = base + strrep(r, " ", "%20"); %#ok<AGROW>
    end
    assert(~isempty(name), 'shx:listICGEM:noFiles', ...
        ['No files parsed from %s - check the Series path against ' ...
         'the temporal catalogue (T.path).'], src);
    T = table(name, url);
    info.note = "Single files download directly (websave on T.url).";
else
    % full temporal catalogue: every series page linked as /sp/<path>.
    % NOTE: no nested capture groups - MATLAB's regexp returns tokens
    % for the OUTERMOST parentheses only (documented behavior, unlike
    % PCRE; caused an index crash in v2.4.1). One flat group.
    m = regexp(html, 'href="/sp/([^"]+)"', 'tokens');
    grp = strings(0, 1); ctr = strings(0, 1); ser = strings(0, 1);
    pth = strings(0, 1); url = strings(0, 1); zp = strings(0, 1);
    seen = strings(0, 1);
    for k = 1:numel(m)
        pk = strrep(string(m{k}{1}), "&amp;", "&");
        if any(seen == pk), continue; end
        seen(end+1, 1) = pk; %#ok<AGROW>
        segs = split(pk, "/");
        segs = segs(strlength(segs) > 0);
        if numel(segs) < 2, continue; end
        grp(end+1, 1) = segs(1); %#ok<AGROW>
        if numel(segs) >= 3
            ctr(end+1, 1) = segs(2); %#ok<AGROW>
            ser(end+1, 1) = join(segs(3:end), "/"); %#ok<AGROW>
        else
            % two-level groups (02_COST-G_, 04_SLR_): center from group
            ctr(end+1, 1) = regexprep(segs(1), '^[\d_]+|[\d_]+$', ''); %#ok<AGROW>
            ser(end+1, 1) = segs(2); %#ok<AGROW>
        end
        pth(end+1, 1) = pk; %#ok<AGROW>
        enc = strrep(pk, " ", "%20");
        url(end+1, 1) = base + "/sp/" + enc; %#ok<AGROW>
        zp(end+1, 1) = base + "/getseries/" + enc; %#ok<AGROW>
    end
    assert(~isempty(grp), 'shx:listICGEM:parseFailed', ...
        'No series parsed from %s - page layout may have changed.', src);
    T = table(grp, ctr, ser, pth, url, zp, 'VariableNames', ...
        {'group', 'center', 'series', 'path', 'url', 'zip'});
    info.note = ['Full temporal catalogue from the series-page links. ' ...
        'zip = whole-series download (hundreds of MB); pass a path ' ...
        'into Series= to list single files. For ITSG monthlies ' ...
        'shx.fetchITSG remains the convenient route.'];
end
end
