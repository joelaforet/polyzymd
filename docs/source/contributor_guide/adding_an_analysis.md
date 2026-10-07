# Add an analysis

An analysis in PolyzyMD is a Python function of an MDAnalysis `AtomGroup` or
`Universe`. The study API runs it on every replicate of every condition:
`Study.timeseries` calls it on each production frame, and
`Study.per_replicate` calls it once per replicate. PolyzyMD stores each
replicate's result with a record of the code, arguments and input files that
produced it, and computes the intervals and tests between conditions with the
replicate as the sampling unit. There is no base class to subclass and no
registry to edit.

{doc}`../how_to/study_api` shows how to load a study, measure on every frame,
compare with a reference structure, turn a replicate into one value, compare
conditions and draw figures. {doc}`../reference/study_api` gives every
signature.

:::{admonition} Environment Setup
:class: tip

Run the examples and tests below in the PolyzyMD analysis pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## Write the function

```python
import polyzymd as pz


def end_to_end(chain):
    """Distance in angstrom between the first and last atom of an AtomGroup."""
    import numpy as np

    return float(np.linalg.norm(chain.positions[-1] - chain.positions[0]))


study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA": "sbma/config.yaml"},
    equilibration="10ns",
)
series = study.timeseries(end_to_end, pz.select("chainid C"), unit="A")
print(series.reduce("mean").compare().to_agent_text())
```

- Arguments given as `pz.select("...")` become the `AtomGroup` of that
  selection in each replicate's `Universe`; `pz.universe()` passes the
  `Universe` itself.
- Import heavy dependencies inside the function, as above: a function in
  `polyzymd.analyses.functions` must not slow down `import polyzymd`.
- A `Study.timeseries` function returns one number for the current frame. A
  `Study.per_replicate` function also receives the keyword `frames`, the
  replicate's production frame indices, and returns one number, or a NumPy
  array with one entry per label (for example per residue) when the call
  gives `labels=`.

## Ship it with PolyzyMD

Add the function to `src/polyzymd/analyses/functions.py` when it is general
enough for other studies. Every function there is written against the same
interface, so the existing ones (`radius_of_gyration`, `rmsf`,
`hydrogen_bonds`, `residue_occlusion`, ...) are working examples. List it in
{doc}`../reference/analysis_functions`.

To make it a `polyzymd analyze NAME` analysis as well, add `NAME` with its
settings and their defaults to `FUNCTION_ANALYSES` in
`src/polyzymd/analyses/protocols.py`, and have `_analyze_function` (or an
`_analyze_<name>` function it calls) run it through the study and return the
`ProtocolReport`. `polyzymd analyze` refuses any name that is not in
`FUNCTION_ANALYSES`, and `--submit` runs any name that is.

## Test it

- Test the function on a small `MDAnalysis.Universe` built in the test, and
  compare it with a value computed by hand or with the MDAnalysis analysis it
  wraps.
- Test the `polyzymd analyze` path with a real config and a few frames; the
  tests in `tests/analyses/` (for example `test_rmsf.py` and
  `test_hydrogen_bonds_analyze.py`) show the fixtures.
- Add a verification page under `docs/source/explanation/` when the analysis
  reproduces a published or independently computed result.

## See also

- {doc}`../reference/study_api`: the study API
- {doc}`../reference/analysis_functions`: the shipped functions
- {doc}`../reference/analysis_protocol_report`: the fields of a `ProtocolReport`
