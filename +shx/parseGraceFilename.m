function meta = parseGraceFilename(filename)
%PARSEGRACEFILENAME Product type and mid-epoch from a GRACE L2 filename.
%
%   META = shx.parseGraceFilename(FILENAME) parses the standard GRACE /
%   GRACE-FO Level-2 naming convention
%       PPP-2_YYYYDDD-YYYYDDD_*.gfc[.gz]
%   e.g. GSM-2_2024032-2024060_GRFO_UTCSR_BA01_0600.gfc, where PPP is one
%   of GSM, GAA, GAB, GAC, GAD and YYYYDDD are year + day-of-year of the
%   coverage span. Best-effort: fields are NaN/"unknown" when the pattern
%   does not match (static models, renamed files); callers can always
%   override via explicit options.
%
%   Inputs
%     filename  char/string   path or bare filename
%   Outputs (struct META)
%     productType  string   "GSM"|"GAA"|"GAB"|"GAC"|"GAD"|"unknown"
%     epoch        double   mid-epoch, decimal years (NaN if unparsed)
%     epochStart   double   coverage start, decimal years (NaN)
%     epochStop    double   coverage end, decimal years (NaN)
%
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     meta       struct: center, product, release, epochStart/Stop (decimal years), raw tokens
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    filename {mustBeTextScalar}
end

[~, base, ext] = fileparts(char(filename));
if strcmpi(ext, '.gz'), [~, base] = fileparts(base); end

meta.productType = "unknown";
meta.epoch = NaN; meta.epochStart = NaN; meta.epochStop = NaN;

tok = regexp(base, '^(G[SA][MABCD])-2_(\d{4})(\d{3})-(\d{4})(\d{3})', ...
    'tokens', 'once');
if isempty(tok)
    % product prefix alone (e.g. renamed monthly file keeping the PPP tag)
    tok2 = regexp(base, '^(GSM|GAA|GAB|GAC|GAD)\>', 'tokens', 'once');
    if ~isempty(tok2), meta.productType = string(tok2{1}); end
    % ITSG daily Kalman names: ..._YYYY-MM-DD (v2.4.1)
    tokD = regexp(base, '_(\d{4})-(\d{2})-(\d{2})$', 'tokens', 'once');
    if ~isempty(tokD)
        y = str2double(tokD{1}); mo = str2double(tokD{2});
        dd = str2double(tokD{3});
        if mo >= 1 && mo <= 12 && dd >= 1 && dd <= 31
            doy = datenum(y, mo, dd) - datenum(y, 1, 1) + 1; %#ok<DATNM>
            meta.epochStart = decYear(y, doy) - 0.5/365.25;
            meta.epochStop  = decYear(y, doy) + 0.5/365.25;
            meta.epoch      = decYear(y, doy);
            if startsWith(base, "ITSG"), meta.productType = "GSM"; end
            return
        end
    end
    % ITSG-style monthly names: ..._YYYY-MM (v2.4.1; covers
    % ITSG-Grace2018 and ITSG-Grace_operational). GSM-equivalent level.
    tok3 = regexp(base, '_(\d{4})-(\d{2})$', 'tokens', 'once');
    if ~isempty(tok3)
        y = str2double(tok3{1}); mo = str2double(tok3{2});
        if mo >= 1 && mo <= 12
            meta.epochStart = y + (mo - 1) / 12;
            meta.epochStop  = y + mo / 12;
            meta.epoch      = y + (mo - 0.5) / 12;
            if startsWith(base, "ITSG"), meta.productType = "GSM"; end
        end
    end
    return
end

meta.productType = string(tok{1});
y1 = str2double(tok{2}); d1 = str2double(tok{3});
y2 = str2double(tok{4}); d2 = str2double(tok{5});
meta.epochStart = decYear(y1, d1);
meta.epochStop  = decYear(y2, d2);
meta.epoch = (meta.epochStart + meta.epochStop) / 2;
end

function y = decYear(year, doy)
nd = 365 + double(mod(year,4)==0 && (mod(year,100)~=0 || mod(year,400)==0));
y = year + (doy - 0.5) / nd;
end
