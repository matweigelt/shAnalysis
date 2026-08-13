function mm = normalizeMonthList(months, errId)
%NORMALIZEMONTHLIST "YYYY-MM" strings or numeric years -> unique list.
%   Shared private helper of the fetch family (fetchSINEX,
%   fetchITSGBackground); mirrors the fetchITSG validation exactly:
%   numeric years must lie in 2002..2100, strings must match
%   ^\d{4}-\d{2}$. ERRID is the caller's error-identifier prefix.
%
%   Inputs
%     months  (1 x k string | numeric)  months or years
%     errId   (1 x 1) string  e.g. "shLowLevel:fetchSINEX"
%
%   Outputs
%     mm  (n x 1) string  unique "YYYY-MM" list
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.12.0).
if isnumeric(months)
    mm = strings(0, 1);
    for y = months(:)'
        assert(y >= 2002 && y <= 2100, char(errId + ":badMonth"), ...
            'Year %g outside 2002..2100.', y);
        mm = [mm; compose("%04d-%02d", y, (1:12)')]; %#ok<AGROW>
    end
else
    mm = string(months(:));
    ok = ~cellfun('isempty', regexp(cellstr(mm), '^\d{4}-\d{2}$', 'once'));
    assert(all(ok), char(errId + ":badMonth"), ...
        'Months must be numeric years or "YYYY-MM" strings.');
end
mm = unique(mm);
end
