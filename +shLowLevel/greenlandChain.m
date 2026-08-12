function [gt, rep] = greenlandChain(folder, gravisFolder, opts)
%GREENLANDCHAIN The validated GravIS Greenland ice-trend chain, one call.
%
%   [GT, REP] = shLowLevel.greenlandChain(FOLDER, GRAVISFOLDER, kn=KN)
%   reproduces guide chapters V1-V6 end to end: gravisL2B corrections
%   with the ICE-6G_D GIA rate, EWH synthesis on a global grid, a
%   pixel-wise trend + annual + semiannual fit, and the regularised
%   leakage inversion stopped on the discrepancy principle with the
%   TREND-matched noise level sigma_monthly/sqrt(Sxx) (the V4c lesson).
%   With the tested defaults and the COST-G RL02.1 series this returns
%   -225.7 Gt/yr against the published -231.1 (2.4%; span
%   2002-04..2023-02). Every input is exchangeable for sensitivity
%   studies; every step lands in REP.
%
%   Inputs
%     folder        (1,1) string  monthly GSM-2_*.gfc folder
%     gravisFolder  (1,1) string  GravIS aux folder (see gravisL2B); must
%                   also hold BasinFile unless given as a full path
%
%   Options
%     kn            ([])  load Love numbers (degree,kn) or column
%                   vector - REQUIRED, the empty default errors: no
%                   frame is assumed for you
%     BasinFile     ("basins_GIS.json")  GravIS GeoJSON of the target
%                   basins ([lon lat] rings; converted internally - see
%                   https://gravis.gfz.de/basins/GIS)
%     NeighbourBoxes ([60 84 232 300; 63 67 335 347; 76 81 10 34])
%                   the guide-V2 Canadian Arctic, Iceland and Svalbard
%                   boxes, [latMin latMax lonMin lonMax] rows, lon in
%                   [0, 360)
%                   every region that can hold mass belongs in the
%                   leakage mask; [] disables the union
%     GIAFile       ("GRAVIS-2B_COSTG_0200_GIA_ICE-6G_D_VM5a_0001.gz")
%     Filter        ("gauss445")  filter DECLARED to the inversion and
%                   applied to the synthesis ("none" | "gaussN")
%     SpanEnd       (2023.099)    published-trend span (guide V3)
%     GridStep      (1)           synthesis grid step [deg]
%     NoiseLevel    (NaN)         NaN: derive sigma_trend from the
%                   open-ocean residual RMS and the regressor power
%                   (V4c); a number overrides it
%     OceanMask     ([])          USER-SUPPLIED open-ocean mask for the
%                   sigma_trend policy (logical grid or @(lat,lon)
%                   handle; see oceanRMS - no coastline is assumed).
%                   REQUIRED unless NoiseLevel is numeric. The mask
%                   must be FALSE over land: a pure latitude band lets
%                   continental hydrology inflate sigma and the
%                   inversion stops early (verified: -220.7 instead of
%                   -225.6 Gt/yr). Guide V4c ships a crude
%                   continent-box mask that reproduces the headline
%     MaxIter       (400)
%     TrendGrid     ([])          precomputed struct from a previous
%                   REP (fields trendGrid, lat, lon, sigTrend) - skips
%                   the synthesis/fit stage for sensitivity studies
%     GM            (3.986004415e14), R (6378136.3)
%     Quiet         (false)
%
%   Outputs
%     gt   (1,1) double  ice-mass trend over the basin mask [Gt/yr]
%     rep  (1,1) struct  trendGrid, lat, lon, sigTrend, mask, unionMask,
%          iterations, stoppedBy, m (inverted map), steps, l2b (gravisL2B
%          report), version, created
%
%   Example
%     kn = readmatrix("loadLoveNumbers_Gegout97.txt", FileType="text", ...
%         NumHeaderLines=2);
%     [gt, rep] = shLowLevel.greenlandChain("E:/series/COSTG", ...
%         "E:/GravIS", kn = kn, OceanMask = @(la,lo) abs(la) < 60);
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.8.8).
arguments
    folder (1,1) string
    gravisFolder (1,1) string
    opts.kn double = []
    opts.BasinFile (1,1) string = "basins_GIS.json"
    opts.NeighbourBoxes double = [60 84 232 300; 63 67 335 347; 76 81 10 34]
    opts.GIAFile (1,1) string = "GRAVIS-2B_COSTG_0200_GIA_ICE-6G_D_VM5a_0001.gz"
    opts.Filter (1,1) string = "gauss445"
    opts.SpanEnd (1,1) double = 2023.099
    opts.GridStep (1,1) double {mustBePositive} = 1
    opts.NoiseLevel (1,1) double = NaN
    opts.OceanMask = []
    opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 400
    opts.TrendGrid = []
    opts.GM (1,1) double = 3.986004415e14
    opts.R (1,1) double = 6378136.3
    opts.Quiet (1,1) logical = false
end
[gt, rep] = shLowLevel.gravisRegionChain("greenland", folder, gravisFolder, opts);
end
