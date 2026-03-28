# Architecture Rules

## Module Structure

```
src/polyzymd/
├── cli/          # Click CLI entry point and command groups
├── config/       # Pydantic v2 configuration models
├── builders/     # System construction pipeline
├── simulation/   # OpenMM simulation execution
├── workflow/     # Orchestration layer
├── core/         # Shared base classes and types
├── analyses/     # ★ Plugin system — unified analysis lifecycle
│   ├── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
│   ├── _*.py     #   Private compute layer (calculators, result models)
│   └── *.py      #   Public plugin files (one per analysis type)
├── compare/      # Cross-condition statistics, formatters, config, IO
├── exporters/    # Format converters (GROMACS, etc.)
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── configs/      # Default YAML configuration files
```

### Inside `analyses/`

| Layer | Files | Role |
|-------|-------|------|
| **Plugins** (public) | `rmsf.py`, `contacts.py`, `rg.py`, etc. | One class per analysis type — the extension point |
| **Compute** (private) | `_calculator_*.py`, `_results_*.py` | Per-condition calculators and result models |
| **Shared utilities** | `shared/loader.py`, `shared/alignment.py`, etc. | `TrajectoryLoader`, alignment, statistics |
| **Framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py` | Plugin ABC, auto-discovery, lifecycle runner |

New analysis types are added as plugins in `analyses/`. The private
`_calculator_*.py` modules provide underlying computation that some plugins
delegate to; new plugins can compute directly in `compute_replicate()`.

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
analyzer = ContactAnalyzer.from_config(config)

# From YAML file
config = SimulationConfig.from_yaml("config.yaml")
```

### ABC + Strategy Pattern
Analysis criteria and molecular selectors use abstract base classes:

```python
# Abstract base
class ContactCriteria(ABC):
    @abstractmethod
    def is_contact(self, distance: float) -> bool: ...

# Concrete strategy
class DistanceCutoffCriteria(ContactCriteria):
    def __init__(self, cutoff: float = 4.5):
        self.cutoff = cutoff

    def is_contact(self, distance: float) -> bool:
        return distance <= self.cutoff
```

### Plugin Discovery Pattern
Analysis plugins in `analyses/` are auto-discovered via `pkgutil`:

```python
from polyzymd.analyses import get_analysis, list_analyses

# List all discovered plugins
for name, cls in list_analyses().items():
    print(f"{name}: {cls.__doc__.splitlines()[0]}")

# Get a specific plugin
RMSFAnalysis = get_analysis("rmsf")
```

No registries, no decorators, no explicit imports needed — just drop a `.py`
file in `analyses/` that subclasses `Analysis`.

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

1. Create `src/polyzymd/analyses/<name>.py`
2. Subclass `Analysis` with `name`, `Settings`, `compute_replicate()`, `aggregate()`
3. For default comparison: implement `extract_metrics()` and set `AggregatedResultClass`
4. Optionally implement `plot()`, `format()`
5. Test with `pixi run -e build pytest tests/test_<name>_plugin.py -v`
6. The CLI automatically discovers it via `polyzymd compare run <name>`

See `analyses/base.py` for the full contract, `analysis-module.md` for
detailed patterns, and `docs/source/tutorials/extending_analyses.md` for the
contributor tutorial.

### Adding a new per-condition calculator

1. Create a private `_calculator_<name>.py` module in `analyses/`
2. Define result models inheriting from `BaseAnalysisResult` (in `_results_<name>.py`)
3. Implement the calculator with `from_config()` factory method
4. Create an `analyses/<name>.py` plugin to expose it through the CLI

### Adding comparison statistics or formatters

The `compare/` package provides shared statistical infrastructure. New plugins
should NOT create files in `compare/` — keep plotting and formatting inline
in the plugin's `plot()` and `format()` methods. The `compare/results/`
directory is used by existing plugins for historical result models only.
