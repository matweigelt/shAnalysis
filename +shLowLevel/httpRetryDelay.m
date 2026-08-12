function w = httpRetryDelay(attempt, retryAfter, opts)
%HTTPRETRYDELAY Wait time for a retry: server hint first, else backoff.
%
%   W = shLowLevel.httpRetryDelay(ATTEMPT, RETRYAFTER) returns the wait
%   [s] before retry number ATTEMPT: when the server sent a finite
%   Retry-After value that value wins (capped at MaxDelay); otherwise
%   exponential backoff BaseDelay * 2^(ATTEMPT-1) with up to 25%
%   uniform jitter, capped at MaxDelay. Pure function - the testable
%   half of shLowLevel.httpFetch.
%
%   Inputs
%     attempt    (1 x 1) double  1-based attempt counter
%     retryAfter (1 x 1) double  server Retry-After [s], NaN if absent
%
%   Options
%     BaseDelay (2)  (1 x 1) first backoff delay [s]
%     MaxDelay (60)  (1 x 1) cap on the returned wait [s]
%     Jitter (0.25)  (1 x 1) relative jitter range, 0 disables
%
%   Outputs
%     w  (1 x 1) double  wait before this retry [s]
%
%   Example
%     w = shLowLevel.httpRetryDelay(3, NaN)        % ~8 s (+jitter)
%     w = shLowLevel.httpRetryDelay(1, 30)         % 30 s (server wins)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.8.9).
arguments
    attempt (1,1) double {mustBeInteger, mustBePositive}
    retryAfter (1,1) double
    opts.BaseDelay (1,1) double {mustBePositive} = 2
    opts.MaxDelay (1,1) double {mustBePositive} = 60
    opts.Jitter (1,1) double {mustBeNonnegative} = 0.25
end
if isfinite(retryAfter) && retryAfter >= 0
    w = min(retryAfter, opts.MaxDelay);
    return
end
w = opts.BaseDelay * 2^(attempt - 1);
w = w * (1 + opts.Jitter * rand());
w = min(w, opts.MaxDelay);
end
