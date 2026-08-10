function writeAnimation(ts, filename, opts)
%WRITEANIMATION Export a monthly series as an MP4 map animation.
%
%   shx.writeAnimation(TS, 'ewh.mp4', quantity="ewh", kn=kn) renders one
%   shx.plotSHMap frame per epoch into an MPEG-4 via VideoWriter (base
%   MATLAB), with a FIXED color scale over all epochs (robust 98th
%   percentile across the whole stack - the single most common animation
%   mistake is a per-frame scale).
%
%   Inputs
%     ts        (1,1) shSeries
%     filename  char/string  .mp4 (Windows/macOS) or .avi (all platforms) target
%   Options
%     quantity ("ewh"), kn, hn        as in synthesis
%     lat (-89:2:89), lon (0:2:358)   grid [deg]
%     FrameRate (4), CLim ([] = robust fixed), Coast (true),
%     Units (""), Projection ("plate")
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    ts (1,1) shSeries
    filename {mustBeTextScalar}
    opts.quantity (1,1) string = "ewh"
    opts.kn double = []
    opts.hn double = []
    opts.lat (1,:) double = -89:2:89
    opts.lon (1,:) double = 0:2:358
    opts.FrameRate (1,1) double = 4
    opts.Sidecar (1,1) logical = true
    opts.CLim double = []
    opts.Coast (1,1) logical = true
    opts.Units (1,1) string = ""
    opts.Projection (1,1) string ...
        {mustBeMember(opts.Projection, ["plate","hammer"])} = "plate"
end
T = ts.nEpochs;
grids = cell(T, 1);
for k = 1:T
    g = ts.at(k);
    args = {opts.lat, opts.lon, 'quantity', char(opts.quantity)};
    if ~isempty(opts.kn), args = [args, {'kn', opts.kn}]; end %#ok<AGROW>
    if ~isempty(opts.hn), args = [args, {'hn', opts.hn}]; end %#ok<AGROW>
    grids{k} = shx.shSynthesis(g.C, g.S, g.GM, g.R, args{:});
end
cl = opts.CLim;
if isempty(cl)
    allv = abs(cat(3, grids{:}));
    a = shx.pctile(allv, 98);
    if a == 0, a = 1; end
    cl = [-a, a];
end
% v2.6.0: profile by extension; MPEG-4 is unavailable in Linux
% VideoWriter, so .avi (Motion JPEG, all platforms) is the portable
% choice and CI-safe format
[~, ~, ext] = fileparts(char(filename));
switch lower(ext)
    case '.mp4'
        profs = VideoWriter.getProfiles;
        if ~any(strcmp({profs.Name}, 'MPEG-4'))
            error('shx:writeAnimation:mp4Unsupported', ...
                ['MPEG-4 is not available on this platform (Linux). ', ...
                 'Use an .avi filename instead.']);
        end
        vw = VideoWriter(char(filename), 'MPEG-4');
    case '.avi'
        vw = VideoWriter(char(filename), 'Motion JPEG AVI');
    otherwise
        error('shx:writeAnimation:badExtension', ...
            'Use an .mp4 or .avi filename (got %s).', ext);
end
vw.FrameRate = opts.FrameRate;
open(vw);
closer = onCleanup(@() close(vw));
fig = figure('Visible', 'off', 'Position', [50 50 900 470]);
figCloser = onCleanup(@() close(fig)); %#ok<NASGU>
ax = axes('Parent', fig);
for k = 1:T
    cla(ax);
    shx.plotSHMap(grids{k}, opts.lat, opts.lon, ax = ax, CLim = cl, ...
        Coast = opts.Coast, Units = opts.Units, ...
        Projection = opts.Projection, ...
        Title = sprintf('%s  %.3f', opts.quantity, ts.epochs(k)));
    drawnow;
    writeVideo(vw, getframe(fig));
end
if opts.Sidecar
    writeSidecar(string(filename), struct('nEpochs', ts.nEpochs, ...
        'quantity', char(opts.quantity), 'frameRate', opts.FrameRate));
end
end

function writeSidecar(target, extra)
% provenance JSON sidecar: <target>.provenance.json (v2.7.0)
v = shx.version();
p = struct('tool', sprintf('%s %s (%s)', v.Name, v.Version, v.Date), ...
    'provenance', v.Provenance, ...
    'created', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'matlab', sprintf('%s / %s', version, computer), ...
    'file', char(target));
fn = fieldnames(extra);
for k = 1:numel(fn), p.(fn{k}) = extra.(fn{k}); end
fid = fopen(char(string(target) + ".provenance.json"), 'w');
fprintf(fid, '%s\n', jsonencode(p, 'PrettyPrint', true));
fclose(fid);
end
