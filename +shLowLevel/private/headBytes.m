function b = headBytes(url, timeoutSec, proxy)
%HEADBYTES Content-Length of a URL via HTTP HEAD (package-private).
%
%   B = headBytes(URL, TIMEOUT, PROXY) returns the announced download
%   size in bytes, or NaN when the server does not disclose it or the
%   request fails - callers print the size only when known. Base
%   MATLAB (matlab.net.http); PROXY as in webFetch.
%
%   Inputs
%     url         (1,1) string  target URL
%     timeoutSec  (1,1) double  request timeout [s]
%     proxy       (1,1) string  "" or "http://host:port"
%
%   Outputs
%     b          (1,1) double  Content-Length [bytes], NaN if unknown
%
%   Example
%     b = headBytes("https://example.com/f.gfc", 10, "");
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-10 (v3.1.0).
b = NaN;
try
    opt = matlab.net.http.HTTPOptions('ConnectTimeout', timeoutSec);
    if strlength(proxy) > 0
        opt.ProxyURI = matlab.net.URI(proxy);
    end
    req = matlab.net.http.RequestMessage('HEAD');
    resp = req.send(matlab.net.URI(url), opt);
    f = resp.getFields('Content-Length');
    if ~isempty(f)
        b = str2double(f(1).Value);
    end
catch
end
end
