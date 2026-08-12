# tools/ci - CI template for embedding shAnalysis in another repository

`shanalysis-ci.yml` is a GitHub-Actions template for projects that vendor
this toolbox under `3rdParty/shAnalysis` (note the `paths:` filters and
the `WORKDIR` env). It is NOT the toolbox's own CI - that lives in
`.github/workflows/ci.yml` and is the maintained reference (job timeout,
gate order); adapt this template from there when embedding.
