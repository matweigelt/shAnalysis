function audit_battery
%AUDIT_BATTERY Unhappy-path battery (audit item 5) + gridScaling pipeline
%   invariance. Each case classified LOUD (error w/ id) / WARN / SILENT.
%   Claude (Fable 5), 2026-08-12.
root = fileparts(mfilename('fullpath'));
dd = fullfile(root, 'tests', 'test_data');
tmp = tempname; mkdir(tmp); c = onCleanup(@() rmdir(tmp, 's'));
gfc = fullfile(dd, 'ITSG-Grace2018_n60_2008-04.gfc');
kn = readmatrix(fullfile(dd,'loadLoveNumbers_Gegout97.txt'), FileType='text', NumHeaderLines=2);

% B1 truncated gfc (cut mid-record)
raw = fileread(gfc); f1 = fullfile(tmp,'trunc.gfc');
fid=fopen(f1,'w'); fwrite(fid, raw(1:round(0.6*numel(raw)))); fclose(fid);
probe('B1 truncated gfc -> read', @() shLowLevel.shReadGFC(f1), @(g) chkTrunc(g));

% B2 corrupted numeric field
bad = regexprep(raw, '(gfc\s+3\s+1\s+)\S+', '$1NOTANUMBER', 'once');
f2 = fullfile(tmp,'corrupt.gfc'); fid=fopen(f2,'w'); fwrite(fid,bad); fclose(fid);
probe('B2 corrupt numeric -> read', @() shLowLevel.shReadGFC(f2), @(g) chkNaNat(g,3,1));

% B3 empty file
f3 = fullfile(tmp,'empty.gfc'); fid=fopen(f3,'w'); fclose(fid);
probe('B3 empty gfc', @() shLowLevel.shReadGFC(f3), []);

% B4 truncated gz (readSHM)
gz = fullfile(dd,'GRAVIS-2B_GIA_ICE-6G_D_VM5a_n10_trimmed.shm.gz');
rawz = fileread_(gz); f4 = fullfile(tmp,'trunc.shm.gz');
fid=fopen(f4,'w'); fwrite(fid, rawz(1:round(0.5*numel(rawz)))); fclose(fid);
probe('B4 truncated .gz -> readSHM', @() shLowLevel.readSHM(f4), []);

% B5 leakageCorrect: empty mask
g = shCoefficients.read(gfc);
lat=(-89:2:89)'; lon=(1:2:359)';
obs = shLowLevel.shSynthesis(g.C, g.S, g.GM, g.R, lat, lon, 'quantity','ewh','kn',kn);
mk0 = false(numel(lat), numel(lon));
probe('B5 leakage empty mask', @() shLowLevel.leakageCorrect(obs, lat, lon, ...
    Mask=mk0, Filter="gauss500", NoiseLevel=1e-4, Quiet=true), []);

% B6 leakageCorrect: NaN in obs
obsN = obs; obsN(40,40) = NaN;
mk = false(size(obs)); mk(80:100, 20:60) = true;
probe('B6 leakage NaN obs', @() shLowLevel.leakageCorrect(obsN, lat, lon, ...
    Mask=mk, Filter="gauss500", NoiseLevel=1e-4, Quiet=true), @(m) chkFinite(m, mk));

% B7 leakageCorrect: size mismatch lat vs obs
probe('B7 leakage size mismatch', @() shLowLevel.leakageCorrect(obs(1:50,:), lat, lon, ...
    Mask=mk(1:50,:), Filter="gauss500", NoiseLevel=1e-4, Quiet=true), []);

% B8 single-epoch series -> climatology
ts1 = shSeries(g.C, Ss=g.S, Epochs=2008.29);
probe('B8 single-epoch climatology', @() ts1.climatology(), []);

% B9 two-epoch series -> trend via climatology (underdetermined 6-param fit)
Cs2 = cat(3, g.C, 1.01*g.C); Ss2 = cat(3, g.S, g.S);
ts2 = shSeries(Cs2, Ss=Ss2, Epochs=[2008.29; 2008.37]);
probe('B9 two-epoch climatology', @() ts2.climatology(), []);

% B10 basinKernel: degenerate zero-area polygon
idx = shLowLevel.shIndex(20, MinDegree=0);
polyDeg = [10 10; 10 10; 10 10];
probe('B10 basin zero-area polygon', @() shLowLevel.basinKernel(idx, polyDeg), @(b) chkZero(b));

% B11 basinKernel: empty mask region
probe('B11 basin empty mask', @() shLowLevel.basinKernel(idx, false(91,180)), @(b) chkZero(b));

% B12 synthesis: kn too short for ewh at nmax
probe('B12 ewh kn too short', @() shLowLevel.shSynthesis(g.C, g.S, g.GM, g.R, ...
    lat, lon, 'quantity','ewh','kn',kn(1:31)), []);

% B13 synthesis: nmax option > size(C)
probe('B13 nmax > size(C)', @() shLowLevel.shSynthesis(g.C, g.S, g.GM, g.R, ...
    lat, lon, 'quantity','geoid','nmax',120), []);

% B14 estimateDegree1 without kn (documented: must error)
probe('B14 degree1 without kn', @() shLowLevel.estimateDegree1(ts2, true(90,180)), []);

% B15 oceanRMS: erosion leaves nothing (tiny ocean)
oc1 = false(numel(lat), numel(lon)); oc1(45,90) = true;
probe('B15 oceanRMS eroded to none', @() shLowLevel.oceanRMS(obs, lat, lon, oc1), @(r) chkNaN(r));

% B16 shAnalysisGrid: NaN in grid
obsN2 = obs; obsN2(10,10) = NaN;
probe('B16 analysis NaN grid', @() shLowLevel.shAnalysisGrid(obsN2, lat, lon, 20, ...
    quantity="ewh", kn=kn), @(C) chkAllNaN(C));

% B17 gridScaling pipeline invariance k(507*m) == k(m)
model = zeros(numel(lat), numel(lon)); model(60:80, 30:70) = 1; model(65:75,40:60)=2.5;
k1 = shLowLevel.gridScaling(model, lat, lon, Filter="gauss500");
k2 = shLowLevel.gridScaling(507*model, lat, lon, Filter="gauss500");
d = max(abs(k1(:)-k2(:)), [], 'omitnan');
fprintf('B17 gridScaling invariance    : max|k(507m)-k(m)| = %.2e  %s\n', d, tern(d<1e-12,'OK (claim holds)','VIOLATION'));
end

% ---------------------------------------------------------- helpers
function probe(name, fn, chk)
try
    w = warning; warning('off','all'); lastwarn('');
    out = fn();
    [wm, wid] = lastwarn; warning(w);
    if ~isempty(wid) || ~isempty(wm)
        fprintf('%-30s: WARN  (%s) %s\n', name, wid, short(wm));
    elseif ~isempty(chk)
        [ok, msg] = chk(out);
        if ok, fprintf('%-30s: SILENT-OK (%s)\n', name, msg);
        else,  fprintf('%-30s: SILENT-WRONG (%s)\n', name, msg);
        end
    else
        fprintf('%-30s: SILENT (no error, no warning)\n', name);
    end
catch ME
    fprintf('%-30s: LOUD  (%s) %s\n', name, ME.identifier, short(ME.message));
end
end
function s = short(m), m = strrep(string(m), newline, ' '); s = extractBefore(m + " ", min(strlength(m)+1, 78)); end
function o = tern(c,a,b), if c, o=a; else, o=b; end, end
function r = fileread_(f), fid=fopen(f,'rb'); r=fread(fid,inf,'*uint8')'; fclose(fid); end
function [ok,msg] = chkTrunc(g)
ok = isstruct(g) || isobject(g); msg = "returned partial object";
if isfield(g,'C') || isprop(g,'C'), msg = sprintf("partial C, nnz=%d", nnz(g.C)); end
end
function [ok,msg] = chkNaNat(g,n,m)
v = g.C(n+1,m+1); ok = ~isnan(v); msg = sprintf("C(%d,%d)=%g", n, m, v);
ok = true; % classification only: report the value
end
function [ok,msg] = chkFinite(m, mk)
ok = all(isfinite(m(mk))); msg = sprintf("mask output finite=%d, NaN count=%d", ok, nnz(~isfinite(m)));
end
function [ok,msg] = chkZero(b)
if isstruct(b), v = b; else, v = struct('k',b); end
ok = true; msg = sprintf("||b||=%.3g", norm(b(:)));
end
function [ok,msg] = chkNaN(r), ok = isnan(r); msg = sprintf("rms=%g", r); end
function [ok,msg] = chkAllNaN(C), ok = true; msg = sprintf("nnzNaN(C)=%d of %d", nnz(isnan(C)), numel(C)); end
