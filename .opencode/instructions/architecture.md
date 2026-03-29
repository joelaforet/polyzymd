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
├── analyses/     # ★ Plugin system — unified analysis lifecycle
│   ├── shared/   #   Reusable utilities (TrajectoryLoader, alignment, statistics, etc.)
│   └── <name>/   #   One package per analysis type (all plugins are packages)
├── compare/      # Shared comparison infrastructure (statistics, config, IO, CLI)
├── exporters/    # Format converters (GROMACS, etc.)
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── configs/      # Default YAML configuration files
```

### Inside `analyses/`

| Layer | Files | Role |
|-------|-------|------|
| **Plugins** (public) | `rmsf/`, `contacts/`, `distances/`, etc. | One class per analysis type — the extension point |
| **Private modules** | `<name>/_plotters.py`, `<name>/_results.py`, etc. | Plotting functions, result models, formatters used internally by plugins |
| **Shared utilities** | `shared/loader.py`, `shared/alignment.py`, etc. | `TrajectoryLoader`, alignment, statistics |
| **Shared compute** | `shared/binding_preference.py`, `shared/surface_exposure.py` | Cross-plugin compute (used by contacts, BFE, polymer_affinity) |
| **Framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py` | Plugin ABC, auto-discovery, lifecycle runner |

New analysis types are added as **packages in `analyses/<name>/`**. All existing
plugins are packages (no single-file plugins exist). New plugins can compute
directly in `compute_replicate()` — no private calculator modules needed.

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

No registries, no decorators, no explicit imports needed — just create a
package in `analyses/<name>/` that subclasses `Analysis`.

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

1. Run `polyzymd new-analysis <name>` to scaffold the plugin, OR create `src/polyzymd/analyses/<name>/` manually
2. Subclass `Analysis` with `name`, `Settings`, `compute_replicate()`, `aggregate()`
3. For default comparison: implement `extract_metrics()`
4. Implement `plot()` using `_build_plot_data()` and shared plotting helpers
5. Optionally implement `format()`
6. Test with `pixi run -e build pytest tests/test_<name>_plugin.py -v`
7. The CLI automatically discovers it via `polyzymd compare run <name>`

See `analyses/base.py` for the full contract, `analysis-module.md` for
detailed patterns, and `docs/source/tutorials/extending_analyses.md` for the
contributor tutorial.

### Adding comparison statistics or formatters

The `compare/` package provides shared statistical infrastructure. New plugins
should NOT create files in `compare/` — keep formatting inline in the plugin's
`format()` method. Established plugins extract plotting into `_plotters.py`
modules within their package; new plugins can start with plotting inline in
`plot()` and extract later as complexity grows. The `compare/results/`
directory is used by existing plugins for historical result models only.
