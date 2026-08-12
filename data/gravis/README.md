# data/gravis - GravIS auxiliary data shipped with shAnalysis

Frozen copies (2026-08-12) of the small GravIS files the validated
chains (`shLowLevel.gravisL2B`, `greenlandChain`, `antarcticaChain`,
`twsChain`) use by default. The correction tables ALWAYS trail the
monthly solutions: for epochs after the freeze, fetch fresh copies from
https://isdc-data.gfz.de/grace/GravIS/COST-G/Level-2B/aux_data/ and
point `gravisFolder` there.

| file | content | source |
|---|---|---|
| GRAVIS-2B_..._GRACE+SLR_LOW_DEGREES_0001.dat | C20/C30/C21/S21 | isdc aux_data |
| GRAVIS-2B_..._GEOCENTER_0001.dat | degree-1 (C10/C11/S11) | isdc aux_data |
| GRAVIS-2B_..._MEAN_..._NFIL_0001.gz | 2002-04..2020-03 mean field | isdc aux_data |
| GRAVIS-2B_..._GIA_ICE-6G_D_VM5a_0001.gz | GIA rate, n=256, epoch 2011 | isdc aux_data |
| basins_GIS.json / basins_AIS.json / basins_rivbas.json | GravIS basin polygons ([lon lat] GeoJSON) | https://gravis.gfz.de/basins/{GIS,AIS,rivbas} |

License: the GravIS Level-2B/Level-3 products are published under
CC BY 4.0 via GFZ Data Services (Dahle & Murboeck 2025,
doi:10.5880/COST-G.GRAVIS_02_L2B; Dahle et al. 2025, ESSD 17, 611-631).
Redistribution here is attribution-preserving; cite the DOIs when
publishing results.
