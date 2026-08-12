# tools/audit - evidence scripts of the 2026-08 independent audit

Executable evidence behind the v3.8.5 audit report, the v3.8.6 fixes and
the guide validation chapters (V1-V8). These are AUDIT ARTIFACTS, not
toolbox API: they carry machine-specific data paths (E:\DATAPOOL\...)
and reproduce the published numbers on the acceptance machine as run.
The productionized, path-free equivalents are shLowLevel.greenlandChain,
antarcticaChain and twsChain.

| script | reproduces |
|---|---|
| audit_gravis.m | guide V1-V6 Greenland chain (-225.7 headline + controls) |
| audit_ais.m | guide V7 Antarctica (-125.7, GIA lever, per-basin table) |
| audit_tws.m | guide V8 TWS (11 basins, DDK variants, GIA finding) |
| audit_senskernel.m | sensitivityKernel 2-16% help claim (F-14 withdrawal) |
| audit_degree1.m | estimateDegree1 physical synthetic world (messy controls) |
| audit_battery.m | unhappy-path battery (B1-B17, LOUD/WARN/SILENT triage) |
| audit_alf_ref.json | 50-dps mpmath legendreALF references (n to 300, poles) |
| audit_gw.txt | Gaussian weight reference vector |
