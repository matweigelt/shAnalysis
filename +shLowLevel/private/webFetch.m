function webFetch(url, dest, timeoutSec, proxy)
%WEBFETCH Download URL to DEST; optional per-call proxy (package-private).
%
%   webFetch(URL, DEST, TIMEOUT, "") uses websave (which honours the
%   proxy configured in MATLAB's Web Preferences, if any). With a
%   nonempty PROXY ("http://host:port", authentication via
%   "http://user:pass@host:port"), the download goes through
%   matlab.net.http with an explicit ProxyURI instead - base MATLAB,
%   no toolbox - streaming to disk via a FileConsumer. Used by all
%   shLowLevel.fetch* functions for institutional networks where per-call
%   proxy control is needed.
%
%   Outputs
%     (none)     writes DEST; errors propagate to the caller, which
%                implements the safe-swap semantics
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-10 (v2.7.0).
if strlength(proxy) == 0
    websave(dest, url, weboptions('Timeout', timeoutSec));
    return
end
httpOpts = matlab.net.http.HTTPOptions( ...
    'UseProxy', true, ...
    'ProxyURI', matlab.net.URI(proxy), ...
    'ConnectTimeout', timeoutSec, ...
    'ResponseTimeout', timeoutSec);
req = matlab.net.http.RequestMessage('GET');
consumer = matlab.net.http.io.FileConsumer(char(dest));
resp = req.send(matlab.net.URI(url), httpOpts, consumer);
if resp.StatusCode ~= matlab.net.http.StatusCode.OK
    if isfile(dest), delete(dest); end
    error('shLowLevel:webFetch:httpError', 'HTTP %s for %s', ...
        char(resp.StatusCode), url);
end
end
