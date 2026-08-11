function rmIfFolder(d)
%RMIFFOLDER Remove a folder if it exists; never throw from a destructor.
%   RMIFFOLDER(D) deletes the folder D and its contents when it is
%   there, and does nothing when it is not.
%
%   Test cleanups run from onCleanup destructors, and a destructor that
%   throws produces a warning MATLAB cannot attach to any test - noise
%   in every CI log, and worse, it can mask the real failure that ended
%   the test early. The usual cause is a folder that was never created:
%   `d = tempname` only invents a NAME, so a test whose code path did
%   not reach the mkdir (or which the framework aborted first) leaves
%   rmdir with nothing to remove.
%
%   Inputs
%     d          (1,1) string | char   folder path; may be absent
%   Outputs
%     none - the folder is removed if present
%
%   Example
%     d = tempname;
%     cl = onCleanup(@() rmIfFolder(d));   %#ok<NASGU>
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.2.1).
if isfolder(d)
    [ok, msg] = rmdir(char(d), 's');
    if ~ok
        % still not fatal: warn so a leaked temp folder is visible
        warning('shAnalysis:test:cleanupFailed', ...
            'could not remove %s: %s', char(d), strtrim(msg));
    end
end
end
