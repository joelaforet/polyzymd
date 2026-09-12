---
name: polyzymd-extend
description: Add a new analysis to PolyzyMD in about 40 lines using the observable contract. Use when an analysis does not exist yet and you need to measure something new from trajectories, when asked to "add an analysis", "extend PolyzyMD", or "write a plugin". Do not copy an existing plugin package, and do not write persistence, statistics or plotting code.
---

# Add an analysis to PolyzyMD

## 1. Scaffold it

```bash
pixi run -e analysis polyzymd new-analysis <name> --style contract
```

This writes `src/polyzymd/analyses/<name>.py` and
`tests/analyses/plugins/test_<name>.py`. Edit those two files and nothing else.

## 2. The contract

```python
class MySettings(BaseModel):
    selection: str = "protein"

class My:
    name: ClassVar[str] = "my"
    Settings: ClassVar[type[BaseModel]] = MySettings
    references: ClassVar[tuple[str, ...]] = ("Author 2020, Journal 1:1, doi:...",)

    def compute(self, universe, frames, settings) -> Sequence[Observable]:
        group = universe.select_atoms(settings.selection)
        values = [measure(group) for _ in iter_frames(universe, frames)]
        return [Observable(name="my", kind="mean_of_timeseries", unit="A", values=values)]

MyAnalysis = contract_analysis(My)   # keep this line; discovery needs it
```

`compute` runs once per replicate and returns raw per-frame numbers. Never
average across replicates, test anything, write a file, or import matplotlib.
The framework does all of that from `kind`.

## 3. Pick the kind

| kind | values are | replicate value | reported as |
|---|---|---|---|
| `mean_of_timeseries` | one number per frame | mean | mean, SEM, 95 percent CI over replicates |
| `fluctuation` | one number per frame | sample standard deviation | same, on the fluctuation |
| `fraction` | 0 or 1 per frame | occupancy | same, on the fraction |
| `distribution` | one number per frame | mean, series kept in the NPZ sidecar | same, plus the shape |
| `profile` | one number per index | the vector | per-index mean and SEM, no pairwise test yet |

Every observable states a `unit`. A `profile` also states an `index`, one entry
per value (residue IDs, bin centres).

## 4. Citations

Put the method paper in the module docstring under a NumPy `References`
heading and repeat it in `references`; the comparison artifact carries it.

## 5. The one test to write

Assert a known answer on a synthetic universe (the Rg of a unit cross is 1.0),
then that `aggregate_observables` over three identical replicates gives a zero
SEM and `n_replicates == 3`. The scaffolded test has both; replace the numbers.

```bash
PYTHONPATH=$PWD/src pixi run -e test pytest tests/analyses/plugins/test_<name>.py -q
pixi run -e analysis polyzymd compare run <name> -f comparison.yaml
```

## 6. Check the science

Run the `livecoms-check` skill before opening a pull request. It fails a number
reported without a unit, an uncertainty that does not say what it is, a test on
frames instead of replicates, and an uncited method. The contract handles the
first three as long as you do not aggregate inside `compute`.
