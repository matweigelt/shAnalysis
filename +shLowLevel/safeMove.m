function n = safeMove(src, dst, opts)
%SAFEMOVE Move a file into place, retrying while a scanner still holds it.
%
%   SAFEMOVE(SRC, DST) renames SRC to DST, overwriting DST, and retries
%   with exponential backoff while the move fails. Every shLowLevel.fetch*
%   function downloads to "<target>.part", verifies the fresh file by
%   parsing it, and only then swaps it into place - this function is that
%   swap.
%
%   The retry exists because the swap can fail for a reason that has
%   nothing to do with the toolbox: on Windows an antivirus or a cloud
%   sync client (OneDrive, Dropbox, corporate DLP agents) opens a
%   just-written file to scan it, holding a lock for a fraction of a
%   second to a few seconds. A plain movefile hits that window, fails,
%   and the caller reports a download error for a file that downloaded
%   perfectly. Retrying costs nothing when there is no lock and turns a
%   spurious failure into a short pause when there is one.
%
%   The result is VERIFIED rather than trusted: movefile can report
%   success on some network and synced filesystems while the destination
%   is not yet visible, so each attempt checks that DST exists and SRC
%   is gone before returning.
%
%   Inputs
%     src        (1,1) string   existing file to move (e.g. "x.gfc.part")
%     dst        (1,1) string   destination path; overwritten if present
%
%   Options
%     Retries (5)  number of RETRIES after the first attempt; 0 gives a
%             single attempt and the behaviour of a plain movefile
%     Pause (0.25)  seconds before the first retry. Each further wait
%             doubles, capped at 4 s, so the default schedule is
%             0.25/0.5/1/2/4 s - about 7.75 s of patience in total,
%             which covers the scanner windows observed in practice
%
%   Outputs
%     n          (1,1) double   number of attempts used (1 = the move
%                succeeded immediately, which is the normal case)
%
%   Errors
%     shLowLevel:safeMove:noSource   SRC does not exist
%     shLowLevel:safeMove:noDestFolder   the folder of DST does not
%                exist - a caller bug, reported immediately instead of
%                spending the retry budget on an impossible move
%     shLowLevel:safeMove:locked     still not movable after all retries;
%                the message carries the OS diagnostic and names the
%                usual cause, so the report is actionable
%
%   Example
%     % download, verify, then swap - the fetcher pattern
%     tmpf = "ITSG-Grace2018_n60_2008-04.gfc.part";
%     shLowLevel.shReadGFC(tmpf);              % verify BEFORE the swap
%     shLowLevel.safeMove(tmpf, "ITSG-Grace2018_n60_2008-04.gfc");
%
%   See also movefile, shLowLevel.fetchICGEM, shLowLevel.fetchITSG.
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.1.4).
arguments
    src (1,1) string
    dst (1,1) string
    opts.Retries (1,1) double {mustBeInteger, mustBeNonnegative} = 5
    opts.Pause (1,1) double {mustBeNonnegative} = 0.25
end
if ~isfile(src)
    error('shLowLevel:safeMove:noSource', ...
        'safeMove: source does not exist: %s', src);
end
% a missing destination folder is a caller bug, not a lock: fail at once
% rather than spending the whole backoff budget on an impossible move,
% and do not blame a scanner for it
dstDir = fileparts(dst);
if strlength(dstDir) > 0 && ~isfolder(dstDir)
    error('shLowLevel:safeMove:noDestFolder', ...
        'safeMove: destination folder does not exist: %s', dstDir);
end
msg = '';
for k = 0:opts.Retries
    n = k + 1;
    [ok, msg] = movefile(char(src), char(dst), 'f');
    % trust the filesystem, not the status flag: synced and network
    % volumes have been observed reporting success too early
    if ok && isfile(dst) && ~isfile(src)
        return
    end
    if k < opts.Retries
        pause(min(opts.Pause * 2^k, 4));
    end
end
error('shLowLevel:safeMove:locked', ...
    ['safeMove: could not move\n    %s\n  to\n    %s\n' ...
     'after %d attempt(s). Last message: %s\n' ...
     'On Windows this is usually an antivirus or cloud-sync scanner ' ...
     'still holding the freshly written file. Raise Retries=/Pause=, ' ...
     'exclude the data folder from on-access scanning, or move the ' ...
     'data folder off the synced drive (shLowLevel.dataFolder).'], ...
    src, dst, opts.Retries + 1, strtrim(msg));
end
