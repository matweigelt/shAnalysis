function [files, info] = fetchITSG(months, opts)
%FETCHITSG Download ITSG monthly OR daily solutions from TU Graz.
%
%   FILES = shLowLevel.fetchITSG(2019:2020) downloads all available monthly
%   GSM solutions for the given years into dataFolder/itsg_series
%   (websave, base MATLAB). Already-present files are skipped unless
%   Update (false)=true, which re-downloads them with a safe swap (the fresh
%   file is parse-verified before it replaces the old one); months
%   that do not exist on the server (mission gap 2017-07..2018-05,
%   intra-mission dropouts) are reported in INFO.missing, not errors.
%
%   FILES = shLowLevel.fetchITSG(["2008-04", "2010-11"]) fetches single months.
%
%   FILES = shLowLevel.fetchITSG("2008-04", Product ("monthly")="daily") fetches the DAILY
%   Kalman-smoother solutions (v2.4.1): one .gfc per day (~83 kB,
%   ICGEM format with formal errors, zero_tide), n40 ONLY - Nmax (NaN)
%   resolves to 40 automatically; requesting another Nmax errors. The
%   months/years spec expands to all days (a full year is ~365 files);
%   target dataFolder/itsg_daily (kept separate from the monthly
%   series so shSeries.fromFolder never mixes samplings).
%
%   Release routing (automatic): epochs before 2017-07 come from
%   ITSG-Grace2018, epochs from 2018-06 on from ITSG-Grace_operational.
%
%   Inputs
%     months   numeric years (each expands to 12 months) OR string
%              array "YYYY-MM"
%   Options
%     Dest (fullfile(shLowLevel.dataFolder(), "itsg_series"))  target folder
%              (set a shared location once via shLowLevel.dataFolder(path))
%     Nmax (96)      60 | 96 (server folder monthly_nXX)
%     Timeout (60)   [s] per file
%     Quiet (false)  suppress per-file progress
%     Proxy ("")  per-call proxy URL, e.g. "http://proxy:8080" (empty: MATLAB Web Preferences)
%     BaseURL ("https://ftp.tugraz.at/pub/ITSG/GRACE")  server base URL, or a local mirror folder for offline use
%     Catalog ([])  row numbers of the shLowLevel.listITSG catalogue to fetch completely
%   Outputs
%     files  (1,:) string   paths of all present files (new + existing)
%     info   struct: fetched, updated, skipped, missing (string arrays), url
%
%   Load the result with ts = shSeries.fromFolder(dest). Source:
%   ftp.tugraz.at/pub/ITSG/GRACE (Mayer-Guerr et al.); cite the ITSG
%   release when publishing.
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    months = ""
    opts.Dest (1,1) string = ""
    opts.Nmax (1,1) double = NaN
    opts.Product (1,1) string ...
        {mustBeMember(opts.Product, ["monthly", "daily"])} = "monthly"
    opts.Timeout (1,1) double = 60
    opts.Proxy (1,1) string = ""
    opts.Update (1,1) logical = false
    opts.Quiet (1,1) logical = false
    opts.Release (1,1) string = ""
    opts.BaseURL (1,1) string = "https://ftp.tugraz.at/pub/ITSG/GRACE"
    opts.Catalog (1,:) double = []
end
daily = opts.Product == "daily";
if isnan(opts.Nmax)
    if daily, opts.Nmax = 40; else, opts.Nmax = 96; end
end
if daily
    assert(opts.Nmax == 40, 'shLowLevel:fetchITSG:badNmax', ...
        'Daily Kalman solutions are released at n40 only.');
else
    assert(any(opts.Nmax == [60, 96, 120]), 'shLowLevel:fetchITSG:badNmax', ...
        'Monthly solutions: Nmax must be 60, 96, or 120.');
end
if isempty(opts.Catalog) && (isstring(months) || ischar(months)) ...
        && strlength(string(months)) == 0 && ~isnumeric(months)
    error('shLowLevel:fetchITSG:noSelection', ...
        'Give months/years, "all", or Catalog= (see shLowLevel.listITSG).');
end
base = opts.BaseURL;
localBase = isfolder(base);
% ---- Catalog= selection (rows of shLowLevel.listITSG) fetches folders fully
if ~isempty(opts.Catalog)
    T = shLowLevel.listITSG(BaseURL = opts.BaseURL, Timeout = opts.Timeout);
    files = strings(1, 0); info = struct('fetched', strings(1, 0), ...
        'updated', strings(1, 0), 'skipped', strings(1, 0), ...
        'missing', strings(1, 0));
    for ii = opts.Catalog(:)'
        assert(ii >= 1 && ii <= height(T), 'shLowLevel:fetchITSG:badCatalog', ...
            'Catalog index %g outside 1..%d (see shLowLevel.listITSG).', ii, height(T));
        row = T(ii, :);
        if row.product == "static"
            warning('shLowLevel:fetchITSG:staticSkipped', ...
                'Catalog %d is a static model folder; use shLowLevel.fetchICGEM.', ii);
            continue
        end
        nmaxArg = row.nmax;
        if row.product == "daily", nmaxArg = 40; end
        [fk, ik] = shLowLevel.fetchITSG("all", Release = row.release, ...
            Product = ternary(row.product == "daily", "daily", "monthly"), ...
            Nmax = nmaxArg, Dest = opts.Dest, BaseURL = opts.BaseURL, ...
            Timeout = opts.Timeout, Proxy = opts.Proxy, ...
            Update = opts.Update, Quiet = opts.Quiet);
        files = [files, fk]; %#ok<AGROW>
        info.fetched = [info.fetched, ik.fetched];
        info.updated = [info.updated, ik.updated];
        info.skipped = [info.skipped, ik.skipped];
        info.missing = [info.missing, ik.missing];
    end
    return
end
% ---- months = "all": enumerate every .gfc in the target folder(s)
if (isstring(months) || ischar(months)) && string(months) == "all"
    rels = opts.Release;
    if strlength(rels) == 0
        rels = ["ITSG-Grace2018", "ITSG-Grace_operational"];
    end
    names = strings(1, 0);
    for rr = rels
        if daily
            error('shLowLevel:fetchITSG:allDaily', ...
                'months = "all" with the daily product would fetch tens of GB; give explicit months/years.');
        end
        folder = sprintf("%s/%s/monthly/monthly_n%d", base, rr, opts.Nmax);
        if localBase
            dd = dir(fullfile(char(folder), '*.gfc'));
            names = [names, string({dd.name})]; %#ok<AGROW>
        else
            html = webread(folder + "/", weboptions('Timeout', opts.Timeout));
            tok = regexp(html, 'href="([^"]+\.gfc)"', 'tokens');
            names = [names, string(cellfun(@(t) t{1}, tok, ...
                'UniformOutput', false))]; %#ok<AGROW>
        end
    end
    names = unique(names, 'stable');
    assert(~isempty(names), 'shLowLevel:fetchITSG:emptyFolder', ...
        'No .gfc files found (release/nmax combination present on server?).');
    mm = extractBetween(names, "_n" + string(opts.Nmax) + "_", ".gfc");
    mm = unique(mm(:)', 'stable');
elseif isnumeric(months)
    mm = strings(1, 0);
    for y = months(:)'
        assert(y >= 2002 && y <= 2100, 'shLowLevel:fetchITSG:badMonth', ...
            'Year %g outside 2002..2100.', y);
        mm = [mm, compose("%04d-%02d", y, (1:12)')']; %#ok<AGROW>
    end
else
    mm = string(months(:)');
    ok = ~cellfun('isempty', regexp(cellstr(mm), '^\d{4}-\d{2}$', 'once'));
    assert(all(ok), 'shLowLevel:fetchITSG:badMonth', ...
        'Months must be numeric years or "YYYY-MM" strings.');
end
dest = opts.Dest;
if strlength(dest) == 0
    sub = 'itsg_series';
    if daily, sub = 'itsg_daily'; end
    dest = string(fullfile(shLowLevel.dataFolder(), sub));             % v2.4.1
end
if ~isfolder(dest), mkdir(dest); end

wo = weboptions('Timeout', opts.Timeout);
files = strings(1, 0); fetched = files; skipped = files; missing = files;
updated = files;
if daily && ~opts.Quiet
    nd = 0;
    for m = mm
        nd = nd + eomday(str2double(extractBefore(m, 5)), ...
            str2double(extractAfter(m, 5)));
    end
    fprintf('  daily product: up to %d files (~%.0f MB)\n', nd, nd*0.083);
end
for m = mm
    y = str2double(extractBefore(m, 5));
    mo = str2double(extractAfter(m, 5));
    frac = y + (mo - 0.5) / 12;
    if strlength(opts.Release) > 0
        rel = opts.Release;                      % explicit: no mixing
    elseif frac < 2017.5
        rel = "ITSG-Grace2018";
    else
        rel = "ITSG-Grace_operational";
    end
    if daily
        days = 1:eomday(y, mo);
    else
        days = 0;                               % sentinel: monthly file
    end
    for dd = days
        if daily
            tag = sprintf("%s-%02d", m, dd);
            fn = sprintf("%s_Kalman_n40_%s.gfc", rel, tag);
            url = sprintf("%s/%s/daily_kalman/daily_n40/%04d/%s", ...
                base, rel, y, fn);
        else
            tag = m;
            fn = sprintf("%s_n%d_%s.gfc", rel, opts.Nmax, m);
            url = sprintf("%s/%s/monthly/monthly_n%d/%s", ...
                base, rel, opts.Nmax, fn);
        end
        fp = fullfile(dest, fn);
        present = isfile(fp);
        if present && ~opts.Update
            skipped(end+1) = string(fp); files(end+1) = string(fp); %#ok<AGROW>
            continue
        end
        tmpf = fp + ".part";
        try
            if localBase
                copyfile(char(url), char(tmpf));
            else
                webFetch(url, tmpf, opts.Timeout, opts.Proxy);
            end
            shLowLevel.shReadGFC(tmpf);                % verify BEFORE swap
            movefile(tmpf, fp, 'f');
            files(end+1) = string(fp); %#ok<AGROW>
            if present
                updated(end+1) = string(fp); %#ok<AGROW>
            else
                fetched(end+1) = string(fp); %#ok<AGROW>
            end
            if ~opts.Quiet && (~daily || dd == 1)
                fprintf('  %s %s%s\n', ...
                    ternary(present, 'updated', 'fetched'), fn, ...
                    ternary(daily, ' (and following days...)', ''));
            end
        catch
            if isfile(tmpf), delete(tmpf); end  % partial download
            if present
                % failed refresh: the existing file stays authoritative
                skipped(end+1) = string(fp); files(end+1) = string(fp); %#ok<AGROW>
            else
                missing(end+1) = tag; %#ok<AGROW>
                if ~opts.Quiet && ~daily
                    fprintf('  missing %s (gap or dropout)\n', m);
                end
            end
        end
    end
end
if daily && ~opts.Quiet && ~isempty(missing)
    fprintf('  %d days missing on the server (gaps/dropouts)\n', ...
        numel(missing));
end
info = struct('fetched', fetched, 'updated', updated, 'skipped', skipped, ...
    'missing', missing, 'dest', dest, 'base', base, 'product', opts.Product);
end

function s = ternary(tf, a, b)
if tf, s = a; else, s = b; end
end
