# Architecture Rules

## Module Structure

```
src/polyzymd/
├── cli/          # Click CLI entry point, command groups, scaffold generator
├── config/       # Pydantic v2 configuration models
├── builders/     # System construction pipeline
├── simulation/   # OpenMM simulation execution
├── workflow/     # Orchestration layer
├── core/         # Shared base classes and types
├── analyses/     # ★ Study API (functions over replicate universes) and the plugin framework being removed
│   └── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
├── exporters/    # Format converters (GROMACS, etc.)
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── templates/    # Example YAML configs and project templates
```

### Inside `analyses/`

| Layer | Files | Role |
|-------|-------|------|
| **Study API** (public) | `study.py`, `timeseries.py`, `functions.py`, `reference.py`, `figures.py`, `protocols.py` | Replicates as MDAnalysis universes, shipped measurement functions, statistics, reports and `polyzymd analyze`; where new measurements go |
| **Shared utilities** | `shared/loader.py`, `shared/window.py`, etc. | `TrajectoryLoader`, frame windows, statistics |
| **Plugin framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py`, `mda/`, `_framework/` | The `Analysis` lifecycle behind `polyzymd compare`; no shipped analysis is a plugin any more, and it is being removed |

Every shipped analysis is a function in `functions.py` run through the study
API: polymer-protein contacts (`residue_occlusion`, `residue_contacts`,
`contact_lifetimes`) by `protocols._analyze_contacts`, hydrogen bonds
(`hydrogen_bonds`, `hbond_lifetimes`, `residue_hbond_occupancy`,
`residue_pair_hbond_occupancy`) by `protocols._analyze_hydrogen_bonds`. Do not
add plugins; `polyzymd.analyses.base` stays the import surface of the framework
until it is removed.

## Chain Convention (Critical)

All systems use a standardized chain-ID mapping:

| Chain | Role | Example |
|-------|------|---------|
| A | Protein (enzyme) | Lipase, protease |
| B | Substrate (small molecule) | p-nitrophenyl butyrate |
| C | Polymer (conjugate) | PEG, polyacrylamide |
| D+ | Solvent, ions, others | Water, Na+, Cl- |

This convention is used throughout selections, analysis, and visualization.
Always use `chainid A`, `chainid B`, etc. in MDAnalysis selections.

## Design Patterns

### Factory Pattern
All major classes support construction from config objects or YAML files:

```python
# From config object
engine = OpenMMEngine.from_config(config)

# From YAML file
config = SimulationConfig.from_yaml("config.yaml")
```

### ABC + Strategy Pattern
Molecular selectors and molecule chargers use abstract base classes:

```python
# Abstract base (analyses/shared/selectors/base.py)
class MolecularSelector(ABC):
    @abstractmethod
    def select(self, universe: Universe) -> SelectionResult: ...

# Concrete strategy (analyses/shared/selectors/polymer.py)
class PolymerChains(MolecularSelector):
    def select(self, universe: Universe) -> SelectionResult: ...
```

### Function Pattern
A measurement is a function of MDAnalysis atom groups; the study supplies the
replicate universes, the production frames, the records and the statistics:

```python
import polyzymd as pz
from polyzymd.analyses.functions import radius_of_gyration

study = pz.Study.from_configs({"A": "A/config.yaml", "B": "B/config.yaml"}, equilibration="10ns")
series = study.timeseries(radius_of_gyration, pz.select("protein"), unit="A", name="rg")
print(series.reduce("mean").compare().to_agent_text())
```

### Config Pattern
Pydantic v2 `BaseModel` subclasses with validators:

```python
class AnalysisConfig(BaseModel):
    cutoff: float = Field(default=4.5, gt=0)
    n_workers: int = Field(default=4, ge=1)

    @model_validator(mode="after")
    def validate_config(self) -> Self:
        # Cross-field validation
        return self
```

## Extension Points

### Adding a new analysis type (primary path)

1. Write a function in `analyses/functions.py`: per frame (atom groups in, one
   number out) for `Study.timeseries`, or per replicate (with `frames`) for
   `Study.per_replicate`. Call MDAnalysis or MDTraj; do not reimplement them.
2. Test it on placed geometry with a known answer, and on the synthetic OpenMM
   run directories of `tests/_support/analysis_testkit.py`.
3. To ship it in `polyzymd analyze`, add a `FUNCTION_ANALYSES` entry and an
   `_analyze_<name>` in `protocols.py`, with its figures, and document it in
   `docs/source/reference/analysis_functions.md` and a how-to page.
4. **Test**: `pixi run -e test pytest tests/analyses -v -k <name>`

See `analysis-module.md` for detailed patterns and
`docs/source/explanation/analysis_api.md` for the study API.

### Adding comparison statistics

Cross-replicate statistics live in `analyses/shared/statistics.py` and
`analyses/shared/inferential_statistics.py`, used by
`ReplicateValues.summary()` and `compare()` in `analyses/timeseries.py`.
