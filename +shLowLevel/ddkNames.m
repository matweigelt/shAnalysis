function names = ddkNames()
%DDKNAMES Released Wbd filenames for DDK1..DDK8 (strong -> weak).
%   NAMES = shLowLevel.ddkNames() - single source of truth for the mapping
%   (verified against the strawpants/GRACE-filter repository listing).
%   Claude (Fable 5), 2026-08-07 (v2.4.1).
%   Outputs
%     names      (1 x 8) string   Wbd file names for DDK1..DDK8
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
names = [ ...
    "Wbd_2-120.a_1d14p_4"   % DDK1
    "Wbd_2-120.a_1d13p_4"   % DDK2
    "Wbd_2-120.a_1d12p_4"   % DDK3 (also shipped in tests/test_data)
    "Wbd_2-120.a_5d11p_4"   % DDK4
    "Wbd_2-120.a_1d11p_4"   % DDK5
    "Wbd_2-120.a_5d10p_4"   % DDK6
    "Wbd_2-120.a_1d10p_4"   % DDK7
    "Wbd_2-120.a_5d9p_4"];  % DDK8
end
