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
├── analysis/     # Per-condition analysis calculators (compute layer)
├── analyses/     # ★ Plugin system — unified analysis lifecycle
├── compare/      # Statistics, formatters, plotters, config, IO
├── exporters/    # Format converters (GROMACS, etc.)
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── configs/      # Default YAML configuration files
```

### `analysis/` vs `analyses/` (important distinction)

| Package | Role |
|---------|------|
| `analysis/` | Per-condition calculators, results, aggregation — the **compute** layer. Contains `RMSFCalculator`, `DistanceCalculator`, `ParallelContactAnalyzer`, etc. |
| `analyses/` | **Plugin system** — wraps `analysis/` calculators into a unified lifecycle (compute → aggregate → compare → plot → format). **Primary extension point for contributors.** |

New analysis types are added as plugins in `analyses/`. The `analysis/` package
provides the underlying computation that plugins delegate to.

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
3. Optionally implement `extract_metrics()`, `compare()`, `plot()`, `format()`
4. Test with `pixi run -e build pytest tests/ -v -k <name>`
5. The CLI automatically discovers it via `polyzymd compare run <name>`

See `analyses/base.py` for the full contract and `analysis-module.md` for
detailed patterns.

### Adding a new per-condition calculator

1. Create a new module under `analysis/` (e.g., `analysis/hbonds/`)
2. Define result models inheriting from `BaseAnalysisResult`
3. Implement the calculator with `from_config()` factory method
4. Create an `analyses/<name>.py` plugin to expose it through the CLI

### Adding comparison statistics or formatters

The `compare/` package provides shared infrastructure:
- `compare/statistics.py` — statistical functions (t-tests, ANOVA, Cohen's d)
- `compare/formatters.py` — CLI output formatters
- `compare/plotters/` — matplotlib visualization
- `compare/config.py` — `ComparisonConfig`, plot settings
