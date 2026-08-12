function info = httpFetch(url, dest, opts)
%HTTPFETCH Download a URL to a file with Retry-After-aware retries.
%
%   INFO = shLowLevel.httpFetch(URL, DEST) downloads URL to the file
%   DEST. On transient failures (HTTP 429, 500, 502, 503, 504, or a
%   network error) it retries with exponential backoff and jitter; when
%   the server sends a Retry-After header (seconds or an HTTP-date),
%   that value takes precedence over the backoff. Built on
%   matlab.net.http because websave discards the response headers on
%   error - which is exactly why polite retrying was never possible
%   with it. Motivated by the 503 bursts the GravIS portal returns on
%   bulk fetches (guide V7/V8 data); all toolbox fetchers route their
%   downloads through this function.
%
%   Inputs
%     url   (1 x 1) string  the URL to download (https recommended)
%     dest  (1 x 1) string  target file path (parent folder must exist)
%
%   Options
%     MaxTries (4)     (1 x 1) total attempts before giving up
%     BaseDelay (2)    (1 x 1) first backoff delay [s]; attempt k waits
%                      BaseDelay * 2^(k-1) plus up to 25% jitter
%     MaxDelay (60)    (1 x 1) cap on any single wait [s], also applied
%                      to server-sent Retry-After values
%     Timeout (60)     (1 x 1) per-attempt timeout [s]
%     Quiet (true)     (1 x 1) suppress the per-retry progress line
%
%   Outputs
%     info  (1 x 1) struct with fields
%       status  (1 x 1) double  final HTTP status code
%       tries   (1 x 1) double  attempts used
%       waited  (1 x 1) double  total seconds slept between attempts
%       bytes   (1 x 1) double  bytes written to DEST
%
%   Error identifiers
%     shLowLevel:httpFetch:failed     all attempts exhausted (message
%                                     carries the last status/error)
%     shLowLevel:httpFetch:badStatus  non-retryable HTTP error (4xx
%                                     other than 429)
%
%   Example
%     info = shLowLevel.httpFetch("https://isdc-data.gfz.de/...", ...
%         "C:/data/aux.dat", MaxTries = 5);
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.8.9).
arguments
    url (1,1) string
    dest (1,1) string
    opts.MaxTries (1,1) double {mustBeInteger, mustBePositive} = 4
    opts.BaseDelay (1,1) double {mustBePositive} = 2
    opts.MaxDelay (1,1) double {mustBePositive} = 60
    opts.Timeout (1,1) double {mustBePositive} = 60
    opts.Quiet (1,1) logical = true
end
import matlab.net.http.RequestMessage
import matlab.net.http.HTTPOptions
waited = 0; lastStatus = NaN; lastErr = "";
for attempt = 1:opts.MaxTries
    try
        req = RequestMessage('GET');
        ho = HTTPOptions(ConnectTimeout = opts.Timeout, ...
            ResponseTimeout = opts.Timeout);
        resp = req.send(matlab.net.URI(url), ho);
        code = double(resp.StatusCode);
        lastStatus = code;
        if code == 200
            body = resp.Body.Data;
            if ischar(body) || isstring(body)
                fid = fopen(dest, 'w'); fwrite(fid, char(body)); fclose(fid);
            else
                fid = fopen(dest, 'wb'); fwrite(fid, body, 'uint8'); fclose(fid);
            end
            d = dir(dest);
            info = struct('status', code, 'tries', attempt, ...
                'waited', waited, 'bytes', d.bytes);
            return
        elseif any(code == [429 500 502 503 504])
            ra = getRetryAfter(resp);
            w = shLowLevel.httpRetryDelay(attempt, ra, ...
                BaseDelay = opts.BaseDelay, MaxDelay = opts.MaxDelay);
            if attempt < opts.MaxTries
                if ~opts.Quiet
                    fprintf('  httpFetch: HTTP %d, retry %d/%d in %.1f s\n', ...
                        code, attempt, opts.MaxTries, w);
                end
                pause(w); waited = waited + w;
            end
        else
            error('shLowLevel:httpFetch:badStatus', ...
                'HTTP %d for %s (not retryable).', code, url);
        end
    catch ME
        if ME.identifier == "shLowLevel:httpFetch:badStatus", rethrow(ME); end
        lastErr = string(ME.message);
        w = shLowLevel.httpRetryDelay(attempt, NaN, ...
            BaseDelay = opts.BaseDelay, MaxDelay = opts.MaxDelay);
        if attempt < opts.MaxTries, pause(w); waited = waited + w; end
    end
end
error('shLowLevel:httpFetch:failed', ...
    'giving up on %s after %d attempts (last status %g%s).', ...
    url, opts.MaxTries, lastStatus, ...
    ternStr(strlength(lastErr) > 0, "; " + lastErr, ""));
end

function ra = getRetryAfter(resp)
ra = NaN;
h = resp.getFields('Retry-After');
if ~isempty(h)
    v = string(h(1).Value);
    n = str2double(v);
    if isfinite(n)
        ra = n;
    else
        try
            ra = max(0, seconds(datetime(v, InputFormat = ...
                'eee, dd MMM yyyy HH:mm:ss zzz', TimeZone = 'UTC') ...
                - datetime('now', TimeZone = 'UTC')));
        catch
        end
    end
end
end

function s = ternStr(c, a, b)
if c, s = a; else, s = b; end
end
