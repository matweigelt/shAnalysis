function [Ct, St] = shEvalGFCT(model, epoch)
%SHEVALGFCT Evaluate time-variable Stokes coefficients at a given epoch.
%
%   [CT, ST] = SHEVALGFCT(MODEL, EPOCH) reconstructs
%       C(t) = C0 + trnd*(t-t0) + sum_k [acos_k*cos(2*pi*(t-t0)) ...
%                                         + asin_k*sin(2*pi*(t-t0))]
%   (and analogously for S) from MODEL.C / MODEL.S (constant part) and
%   MODEL.variableTerms (from SHREADGFC), at EPOCH (same units/convention
%   as the t0 values in the file -- typically decimal years).
%
%   ASSUMPTION: 'acos'/'asin' terms are treated as ANNUAL (period = 1,
%   same units as t0/epoch). This matches the common ICGEM gfct
%   convention, but some models encode additional sub-annual terms
%   (semi-annual etc.) under different tags -- inspect
%   MODEL.variableTerms(:).type and extend the switch below if your
%   specific model uses other tags. Verify against the model's own
%   documentation before using for production results.
%
%   If MODEL.variableTerms is empty (static gfc file), returns MODEL.C,
%   MODEL.S unchanged (EPOCH is ignored, with a warning).
%
%   ICGEM format 2.0 (header 'format icgem2.0') is fully supported (v2.1):
%   piecewise gfct validity intervals, per-term [t0, t1) activity windows
%   and explicit periods on acos/asin lines; epochs in yyyymmdd(.hhmm)
%   are converted to decimal years (shx.icgemDate2Year). For 1.0 files
%   the pre-v2.1 behavior is preserved bit-for-bit when epochs are given
%   as decimal years; yyyymmdd epochs are now converted, and an 8-column
%   acos/asin line is read as EIGEN-style 'period in the last column'.
%
%   Claude (Sonnet 4.6), 2026-07-11; ICGEM 2.0 support:
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     C          (nmax+1 x nmax+1) double   coefficients evaluated at the epoch
%     S          (nmax+1 x nmax+1) double   sine coefficients at the epoch
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

Ct = model.C;
St = model.S;

if isempty(model.variableTerms)
    warning('shEvalGFCT:staticModel', ...
        'MODEL has no time-variable terms (static gfc file) -- returning constant C, S.');
    return;
end

vt = model.variableTerms;
isV2 = any(strcmp({vt.type}, 'gfct'));

if isV2
    % ---- ICGEM 2.0: piecewise-defined coefficients. For each (n,m) the
    % gfct piece whose validity interval [t0, t1) contains EPOCH is the
    % constant part; trnd/acos/asin terms apply only while EPOCH lies in
    % their own interval, referenced to their own t0. Requesting an epoch
    % outside every piece raises shEvalGFCT:epochOutside.
    isP = strcmp({vt.type}, 'gfct');
    pcs = vt(isP);
    keys = unique([[pcs.n]', [pcs.m]'], 'rows');
    for q = 1:size(keys, 1)
        n = keys(q,1); m = keys(q,2);
        sel = pcs([pcs.n] == n & [pcs.m] == m);
        inIv = arrayfun(@(p) epoch >= p.t0 && epoch < p.t1, sel);
        if ~any(inIv)
            error('shEvalGFCT:epochOutside', ...
                ['Epoch %.4f lies outside every gfct validity interval ' ...
                 'for (n=%d, m=%d).'], epoch, n, m);
        end
        p = sel(find(inIv, 1));
        Ct(n+1,m+1) = p.C;
        St(n+1,m+1) = p.S;
    end
end

for k = 1:numel(vt)
    t = vt(k);
    n = t.n; m = t.m;
    if strcmp(t.type, 'gfct'), continue; end
    if isV2 && ~(epoch >= t.t0 && epoch < t.t1)
        continue;                                % term inactive at EPOCH
    end
    dt = epoch - t.t0;
    period = t.period;
    if isnan(period), period = 1.0; end
    switch t.type
        case 'trnd'
            Ct(n+1,m+1) = Ct(n+1,m+1) + t.C * dt;
            St(n+1,m+1) = St(n+1,m+1) + t.S * dt;
        case 'acos'
            Ct(n+1,m+1) = Ct(n+1,m+1) + t.C * cos(2*pi*dt/period);
            St(n+1,m+1) = St(n+1,m+1) + t.S * cos(2*pi*dt/period);
        case 'asin'
            Ct(n+1,m+1) = Ct(n+1,m+1) + t.C * sin(2*pi*dt/period);
            St(n+1,m+1) = St(n+1,m+1) + t.S * sin(2*pi*dt/period);
        otherwise
            warning('shEvalGFCT:unknownType', ...
                'Unrecognized variable term type ''%s'' at n=%d, m=%d -- skipped.', ...
                t.type, n, m);
    end
end

end
