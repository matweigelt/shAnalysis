function [files, st] = fetchFileSet(urls, names, destFolder, opts)
%FETCHFILESET Robust file-set download - the shared fetch-family loop.
%
%   Shared private helper of fetchGAX, fetchITSGSINEX and
%   fetchITSGBackground: skip-if-present, Update= safe re-download,
%   MaxFiles/BudgetSec/MaxFailures caps, a polite inter-download
%   pause, a single capped retry on HTTP 429 (the ICGEM lesson), and
%   the choice between the websave stack and the Retry-After-aware
%   shLowLevel.httpFetch (which itself falls back to websave on
%   transport-level failures).
%
%   Inputs
%     urls, names (n x 1) string  full URLs and local file names
%     destFolder  (1 x 1) string  target folder (created if absent)
%
%   Options (defaults mirror the fetch-family front ends)
%     Update (false), MaxFiles (Inf), BudgetSec (Inf), MaxFailures (5),
%     PauseSec (1.5), RetryAfterCap (30),
%     Downloader ("websave") "websave" | "httpFetch",
%     Quiet (false)
%
%   Outputs
%     files (m x 1) string  local paths (downloaded or already present)
%     st    (1 x 1) struct  nListed, nDownloaded, nSkipped, nFailed,
%           nRemaining (budget/failure cut)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.12.0).
arguments
    urls (:,1) string
    names (:,1) string
    destFolder (1,1) string
    opts.Update (1,1) logical = false
    opts.MaxFiles (1,1) double = Inf
    opts.BudgetSec (1,1) double {mustBePositive} = Inf
    opts.MaxFailures (1,1) double {mustBePositive} = 5
    opts.PauseSec (1,1) double {mustBeNonnegative} = 1.5
    opts.RetryAfterCap (1,1) double {mustBePositive} = 30
    opts.Downloader (1,1) string ...
        {mustBeMember(opts.Downloader, ["websave", "httpFetch"])} = "websave"
    opts.Quiet (1,1) logical = false
end
if ~isfolder(destFolder), mkdir(destFolder); end
t0 = tic;
files = strings(0, 1);
nD = 0; nS = 0; nF = 0; nRem = 0;
for k = 1:numel(urls)
    out = fullfile(char(destFolder), char(names(k)));
    if isfile(out) && ~opts.Update
        nS = nS + 1; files(end+1, 1) = string(out); %#ok<AGROW>
        continue
    end
    if nD >= opts.MaxFiles || toc(t0) > opts.BudgetSec
        nRem = nRem + numel(urls) - k + 1;
        break
    end
    try
        if opts.Downloader == "websave"
            try
                websave(out, urls(k), weboptions('Timeout', 60));
            catch MEi
                if strcmp(MEi.identifier, ...
                        'MATLAB:webservices:HTTP429StatusCodeError')
                    pause(opts.RetryAfterCap);          % once, capped
                    websave(out, urls(k), weboptions('Timeout', 60));
                else
                    rethrow(MEi);
                end
            end
        else
            shLowLevel.httpFetch(urls(k), out);
        end
        nD = nD + 1; files(end+1, 1) = string(out); %#ok<AGROW>
        if ~opts.Quiet && mod(nD, 25) == 0
            fprintf('  %d files\n', nD);
        end
        if opts.PauseSec > 0, pause(opts.PauseSec); end
    catch ME
        nF = nF + 1;
        if ~opts.Quiet
            fprintf('  [fail] %s: %s\n', names(k), ME.identifier);
        end
        if nF >= opts.MaxFailures
            warning('shLowLevel:fetchFileSet:tooManyFailures', ...
                'stopping after %d failures (last: %s); %d files not attempted.', ...
                nF, ME.identifier, numel(urls) - k);
            nRem = nRem + numel(urls) - k;
            break
        end
    end
end
st = struct('nListed', numel(urls), 'nDownloaded', nD, 'nSkipped', nS, ...
    'nFailed', nF, 'nRemaining', nRem);
end
