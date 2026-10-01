# Citing the methods a study used

```{admonition} Status: design
:class: warning
This page is the agreed design of method citations. Nothing on it is
implemented yet. It is planned as one slice after the study-folder slices,
followed by a smaller slice for the simulation side; until then
`polyzymd cite` does not exist.
```

PolyzyMD runs methods and libraries other people wrote: superposition,
solvent-accessible surface area, secondary structure, hydrogen-bond
detection, survival analysis, statistical inefficiency, multiple-comparison
correction. Their authors should get credit whenever a PolyzyMD result is
published. This feature gives a user, at the end of a study, one command
that prints citations, as BibTeX or formatted text, for every method and
library behind every stored result, so they paste one block into a
reference manager and miss no one.

## What the user sees

```bash
polyzymd cite --study lipase_363K                  # every stored result of the study
polyzymd cite --study lipase_363K contacts         # one run
polyzymd cite --study lipase_363K --format text    # a formatted reference list
polyzymd cite results/rg/report.json               # one report, no study needed
```

| Where | What |
|---|---|
| `polyzymd cite` | BibTeX (default), formatted text, or Citation File Format references, deduplicated across runs |
| Reports | `provenance.citations` lists the keys; the agent format adds one line, `cite: <n> methods and libraries; polyzymd cite --study S RUN` |
| `polyzymd study freeze` | Writes `references.bib` into the deposit, adds the methods and libraries to `CITATION.cff` `references`, and lists the keys in `manifest.json` |
| Python | `report.citations`, and `pz.cite(study_or_report, format="bibtex")` |

Citations are read from stored reports, so `polyzymd cite` works on a
downloaded study without its trajectories. One line in each report keeps
the citations out of the way while a study is in progress, when they do not
matter yet.

## What gets cited

| Kind | Cited when | Examples |
|---|---|---|
| Methods | A result was produced by them | Shrake and Rupley (1973) and Tien et al. (2013) for SASA and occlusion contacts; Kabsch and Sander (1983) for DSSP; Best, Hummer and Eaton (2013) for native contacts; Kaplan and Meier (1958) and Royston and Parmar (2013) for contact lifetimes; Theobald (2005) and Liu et al. (2010) for superposition |
| Statistics | They ran on a result | Welch's t test; Benjamini and Hochberg (1995); Cohen's d and Hedges' g; Chodera (2016) for equilibration detection and statistical inefficiency; Grossfield et al. (2018) for the uncertainty conventions |
| Libraries | Their code ran | MDAnalysis (Michaud-Agrawal et al. 2011; Gowers et al. 2016), MDTraj, pymbar, NumPy, SciPy, each as its own documentation asks |
| Upstream registrations | A library that registers citations with duecredit registered them for the code PolyzyMD called | MDAnalysis's own citations, such as Smith et al. (2019) for its hydrogen-bond analysis |
| The study's own functions | They ran | The `cite:` of their `study.yaml` entry, or `@pz.cites(...)` |
| PolyzyMD | Always | `polyzymd.citation`, as the study folders already do |

Libraries are cited the way their own "how to cite" pages ask, which
sometimes means several papers, or a paper per submodule; the bibliography
records which page each entry follows.

## How it works

### One bibliography

`src/polyzymd/citations/references.bib`, shipped with the package, holds
every entry, each with a DOI checked against its publisher. The methods and
references page ({doc}`references`) is generated from it, or checked
against it by a test, so the two cannot disagree. NumPy-style `References`
sections in docstrings stay, for readers of the code.

### Tagging code

- `@cites("shrake1973", "tien2013")` on a function declares the citations of
  what it computes. The declaration is an attribute of the function, so it is
  known without running it.
- `cite("benjamini1995")` inside a code path credits a method only when that
  path runs, such as a correction applied only to comparisons.
- Every shipped function in `analyses/functions.py`, the statistics in
  `analyses/shared/`, the reference and alignment code, and the loader is
  tagged.

### Collecting citations for a result

A collector is open while one analysis runs, and gathers:

1. the declared citations of every function `Study.timeseries` and
   `Study.per_replicate` run, **including when the values are read back from
   storage**, because a cached result was still produced by that function;
2. the citations `cite()` adds while statistics and other code paths run;
3. the upstream registrations (below) for the modules those functions use;
4. the citations of the study's own functions.

The keys go into the report (`provenance.citations`), never into a stored
record, so adding or changing citations never recomputes a result.

### Upstream citations through duecredit

MDAnalysis, like several other packages, registers its citations with
duecredit (`due.cite(Doi(...), path="MDAnalysis.analysis.hydrogenbonds...")`).
These registrations are inert unless the `duecredit` package is installed.
PolyzyMD must still pass them on: when a PolyzyMD method uses an MDAnalysis
module that registers a citation, the user receives that citation too.

- **Primary approach:** PolyzyMD depends on duecredit, a small BSD-licensed
  package, and runs an in-process duecredit collector while analyses run,
  without duecredit's environment variable or its `.duecredit.p` file. After
  the run it reads the registrations whose `path` is a module the tagged
  functions declare they use (`@cites(..., uses=["MDAnalysis.analysis.hydrogenbonds"])`),
  converts them to BibTeX, and adds them to the result.
- **Check:** a test lists every `due.cite` registration in the installed
  MDAnalysis and fails when a module PolyzyMD uses gains or changes a
  citation, so the bibliography is updated with it.
- **Fallback:** if the duecredit API does not allow an in-process collector,
  the registrations of the modules PolyzyMD uses are copied into
  `references.bib`, kept current by the same test.

Which of the two holds is decided at the start of the slice, by trying the
duecredit API on MDAnalysis 2.10.

### The study's own functions

```yaml
analyses:
  lid_opening:
    function: analyses/lid.py:lid_distance
    kind: timeseries
    cite: [doi:10.1021/acs.jctc.5b00043, mcgibbon2015]
```

A key from `references.bib` is used as it is. A DOI not in the bibliography
is resolved to BibTeX from `https://doi.org` by content negotiation; when
the network is not available, a minimal entry holding only the DOI is used,
with a warning saying so. Resolved entries are cached in the study folder,
so they are fetched once.

## Simulation side (a later slice)

`polyzymd build` and `polyzymd run` use OpenMM or GROMACS, the OpenFF
toolkit and Interchange, NAGL or AM1-BCC charges, the protein and
small-molecule force fields, a water model and packmol. Their citations
follow from each condition's resolved config, so `polyzymd study freeze` and
`polyzymd cite --study` add them for every condition without any tagging
of code paths.

## Changes by file

| Part | Where |
|---|---|
| Bibliography, registry, collector, formatting | new `src/polyzymd/citations/` |
| Tags | `analyses/functions.py`, `analyses/shared/` (statistics, inferential statistics, autocorrelation, diagnostics), `analyses/reference.py`, `analyses/shared/loader.py` |
| Collecting around runs, cached or not | `analyses/timeseries.py`, `analyses/protocols.py` |
| `provenance.citations` and the `cite:` line | `analyses/protocols.py`, `reference/analysis_protocol_report.md` |
| `polyzymd cite`, `pz.cite` | new `cli/cite.py`, `polyzymd/__init__.py` |
| `cite:` for the study's own functions, DOI resolution | `analyses/study_file.py`, `analyses/user_functions.py` |
| `references.bib` and `CITATION.cff` references at freeze | `analyses/study_freeze.py`, `analyses/study_metadata.py` |
| Generated or checked references page, how-to, agent files | `docs/source/explanation/references.md`, a new how-to, `AGENTS.md`, the analyze skill |
| Tests | every key used exists and has a DOI; each shipped analysis yields its expected keys; a cached rerun still cites; MDAnalysis's registrations are passed on; `cite` works without trajectories; an unreachable network warns |
