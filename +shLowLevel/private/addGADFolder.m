function [Cs, Ss, nGad] = addGADFolder(Cs, Ss, ep, gadFolder)
%ADDGADFOLDER Add AOD1B GAD monthly means onto a coefficient stack.
%   Shared private helper of oceanChain and obpChain: begin-date-matched
%   on the GAD-2_* filenames, added BEFORE any filtering (the GravIS
%   Level-3 order). Epochs without a matching file are left untouched;
%   the caller reports the coverage.
%
%   Inputs
%     Cs, Ss     (n+1 x n+1 x T) coefficient stack, C(n+1, m+1) indexing
%     ep         (T x 1) epochs [decimal years]
%     gadFolder  (1 x 1) string  folder of GAD-2_*.gfc
%
%   Outputs
%     Cs, Ss  the stack with GAD added where covered
%     nGad    (1 x 1) number of epochs restored
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.11.0).
% add back GAD-2 monthly means, matched on the begin date in the name
df = dir(fullfile(char(gadFolder), 'GAD-2_*.gfc'));
tok = regexp({df.name}, 'GAD-2_(\d{4})(\d{3})-', 'tokens', 'once');
yd = @(y) 365 + double(mod(y,4)==0 & (mod(y,100)~=0 | mod(y,400)==0));
begG = nan(numel(tok), 1);
for k = 1:numel(tok)
    y = str2double(tok{k}{1}); d = str2double(tok{k}{2});
    begG(k) = y + (d-1)/yd(y);
end
nGad = 0; nmax = size(Cs, 1) - 1;
for k = 1:numel(ep)
    [dmin, j] = min(abs(begG - (ep(k) - 15/365)));  % ~mid - half month
    if isempty(dmin) || dmin > 0.05, continue; end
    g = shCoefficients.read(fullfile(df(j).folder, df(j).name));
    nm = min(nmax, size(g.C, 1) - 1);
    Cs(1:nm+1, 1:nm+1, k) = Cs(1:nm+1, 1:nm+1, k) + g.C(1:nm+1, 1:nm+1);
    Ss(1:nm+1, 1:nm+1, k) = Ss(1:nm+1, 1:nm+1, k) + g.S(1:nm+1, 1:nm+1);
    nGad = nGad + 1;
end
end
