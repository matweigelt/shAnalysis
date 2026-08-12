function [gt, rep] = antarcticaChain(folder, gravisFolder, opts)
%ANTARCTICACHAIN The validated GravIS Antarctica ice-trend chain.
%
%   [GT, REP] = shLowLevel.antarcticaChain(FOLDER, GRAVISFOLDER, kn=KN)
%   runs the same chain as greenlandChain with the 25 GravIS AIS
%   drainage-basin polygons as the mask (guide chapter V7). With the
%   tested defaults and COST-G RL02.1 this returns -125.7 Gt/yr against
%   the AWI joint-inversion -146.9 - read V7 before interpreting the
%   14%: it is the documented distance between map-based estimators and
%   the Sasgen forward-modelling class (the official TU Dresden kernel
%   grid shows the same per-basin damping, correlation 0.905), and the
%   GIA correction (+52.2 Gt/yr over the AIS) is not optional here.
%   REP.basins reports the per-basin masses over the 25 polygons.
%
%   Inputs
%     folder        (1,1) string  monthly GSM-2_*.gfc folder
%     gravisFolder  (1,1) string  GravIS aux folder (see gravisL2B)
%
%   Options
%     kn            ([])  load Love numbers (degree,kn) - REQUIRED,
%                   the empty default errors: no frame is assumed
%     BasinFile     ("basins_AIS.json")  GravIS AIS GeoJSON
%                   (https://gravis.gfz.de/basins/AIS; [lon lat] rings,
%                   converted and dateline-recentred internally)
%     NeighbourBoxes ([])         Antarctica has no land neighbours
%     GIAFile       ("GRAVIS-2B_COSTG_0200_GIA_ICE-6G_D_VM5a_0001.gz")
%     Filter        ("gauss445")  as in greenlandChain
%     SpanEnd       (2023.099), GridStep (1)
%     NoiseLevel    (NaN)  sigma_trend policy, see greenlandChain
%     OceanMask     ([])   REQUIRED unless NoiseLevel is numeric; must
%                   be false over land (see greenlandChain)
%     MaxIter       (400), TrendGrid ([])
%     GM            (3.986004415e14), R (6378136.3)
%     Quiet         (false)
%
%   Outputs
%     gt   (1,1) double  total AIS trend over the basin union [Gt/yr]
%     rep  (1,1) struct  as greenlandChain, plus basins (table: name,
%          gt) integrated per GravIS polygon
%
%   Example
%     [gt, rep] = shLowLevel.antarcticaChain("E:/series/COSTG", ...
%         "E:/GravIS", kn = kn, OceanMask = @(la,lo) abs(la) < 60);
%     rep.basins
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.8.8).
arguments
    folder (1,1) string
    gravisFolder (1,1) string
    opts.kn double = []
    opts.BasinFile (1,1) string = "basins_AIS.json"
    opts.NeighbourBoxes double = []
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
[gt, rep] = shLowLevel.gravisRegionChain("antarctica", folder, gravisFolder, opts);
end
