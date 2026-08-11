# shAnalysis developer toolchain

Everything needed to regenerate the documentation and to run the quality
gates. Nothing in here ships to users; the toolbox itself is base MATLAB
only. The scripts are Python 3 and, except where noted, stdlib only.

## Quality gates

All must report **zero findings** before a PR. The first two also run in
CI (`.github/workflows/ci.yml`), before MATLAB starts, inside the
required job.

| Gate | Runs from | Checks |
|---|---|---|
| `tools/help_audit.py` | repo root | every public function and method documents each input, each name-value option **with its default**, an Outputs section with size/type tokens, and an example |
| `tools/doc_sync_audit.py` | repo root | the three documentation sources agree — see below |
| `tools/dev/mlint_lite.py FILE...` | anywhere | MATLAB lint without MATLAB |
| `tools/dev/check_api.py` | needs `api_data.py` | guide snippets against the parsed API |
| `tools/dev/api_audit.py` | needs `api_data.py` | every registered example matches its contract |
| `tools/dev/attribution_sweep.py` | repo root | the provenance stamp on every deliverable |

### What `doc_sync_audit` checks

shAnalysis documents itself three times: the in-file help, `html/`, and
the workflow guide. They drift silently — before v3.1.2 all the other
gates were green while the API reference listed 5 of 12 `fetchICGEM`
options, eleven help pages were unreachable from the Help browser, the
tagged release reported the previous version, and the guide advertised a
call that threw. The gate covers, in order of how badly each bites:

- **SNIPPET** — every call in an `html/*.html` `<pre>` block or a guide
  `code("""...""")` block must match the parsed contract: the function
  exists, the positional count is in range, option names are declared,
  and option **values** are of the declared type. The value check matters:
  `standardChain(TN14="TN-14.txt")` names a real option and still throws,
  because `TN14` is a logical and the path belongs in `TN14File`.
- **APIREF** — `html/apiReference.html` contains every entity and every
  one of its options, i.e. it was regenerated after the last API change.
- **TOC** — `html/helptoc.xml` reaches every page, and points at none
  that do not exist.
- **VERSION** — `Contents.m` (the single source of truth) agrees with
  `CITATION.cff`, the top CHANGELOG section and the API reference title;
  no `Unreleased` section survives a release; line 1 of `Contents.m` is
  still a short product name, because `ver()` prints it verbatim.
- **HELP** — no option is "documented" as `see arguments block`, and
  every default in the house style `Name (default)  description` matches
  the arguments block exactly. Sentinel defaults (`NaN`, `[]`, `""`,
  `table()`, ...) may instead be documented by their effect, e.g.
  `T0 (mean(tYears))`.
- **COVERAGE** — every public entity is mentioned in at least one
  narrative source, not only in the generated reference.

Narrative pages are deliberately **not** required to list every option
of every function; that is what the generated API reference is for. The
gate enforces a coverage floor, not completeness.

## Regenerating the documentation

The generators expect a build copy at `/home/claude/shx_build/shAnalysis`
and write `api_data.py` next to themselves. Order matters:

```bash
cp -r <repo> /home/claude/shx_build/shAnalysis
python3 tools/dev/api_extract.py     # source  -> /home/claude/api_data.py
#   merge docs/apiExamples.json into api_data.py (examples + run flags)
python3 tools/dev/make_apiref.py     # api_data -> html/apiReference.html
python3 tools/dev/make_guide.py      # api_data -> docs/..._guide.pdf
```

- **Never pipe these through `head`.** SIGPIPE kills the script before it
  writes its output file, and the failure is silent — it has shipped a
  reference stripped of all 154 examples once.
- `VERSION` is always parsed from `Contents.m`. Six separate bugs came
  from hardcoded version strings; do not add a seventh.
- `SHX_GUIDE_OUT` overrides the guide's output path.

## Guide figures

`tools/dev/guide_assets/` holds the ten figures the guide embeds. They
live in the repository on purpose: they were previously read from a
container path that did not survive a session, which meant the shipped
PDF could not be rebuilt by anyone. `make_figs.py` regenerates three of
them (`d01`, `d04_gains`, `d04_diff`) from the Python validation port;
the rest were recovered from the shipped PDF.

If you ever recover figures from a PDF again: **the image order within a
page does not match the order they appear in the story.** Identify each
by content and check it against its caption in a rendered rebuild.

## Validation scripts

`validate_*.py` hold the Python (numpy/scipy) reference implementations
behind the pinned MATLAB test values. Every numerical method is validated
here **before** the MATLAB implementation — this is what caught the
column-index bug in the `gfct` variable-term parser.

Developed by Matthias Weigelt with the help of Claude.
