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
│   └── <name>/   #   Analysis plugins (single-file simple modules or packages)
├── exporters/    # Format converters (GROMACS, etc.)
├── data/         # Bundled data files (force fields, templates)
├── utils/        # Shared utilities
└── templates/    # Example YAML configs and project templates
```

### Inside `analyses/`

| Layer | Files | Role |
|-------|-------|------|
| **Plugins** (public) | `rmsf/`, `contacts/`, `rg.py`, etc. | One settings model plus one `compute()` per analysis, the extension point |
| **Private modules** | `_framework/` | Internal framework implementation, not a contributor import target |
| **Shared utilities** | `shared/loader.py`, `shared/alignment.py`, etc. | `TrajectoryLoader`, alignment, statistics |
| **Framework** | `contract.py`, `base.py`, `contract_plots.py`, `discovery.py`, `orchestrator.py`, `mda/` | The contract, the lifecycle, auto-discovery, artifact storage |

An analysis is one module under `analyses/`. It declares a settings model and a
`compute(universe, frames, settings)` that returns `Observable` objects, and it
ends with `contract_analysis()`, which builds the `Analysis` subclass. PolyzyMD
owns `ArtifactStore`, `ConditionArtifact`, `ComparisonArtifact`, ensemble
orchestration, statistics and plotting; a plugin writes none of that.

A plugin imports `Observable`, `iter_frames` and `contract_analysis` from
`polyzymd.analyses.contract`, and anything else from
`polyzymd.analyses.shared`. `polyzymd.analyses.base` holds `Analysis` and the
four framework contexts a plugin receives. Do not import
`polyzymd.analyses._framework` from a plugin.

`polyzymd.analyses.contacts` follows the same facade pattern. The public
`ContactsAnalysis` class remains in `contacts/__init__.py`; artifact handling,
condition filtering, custom comparison, plotting orchestration, result models,
and MDAnalysis helpers live in private `contacts/_*.py` modules.

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

No registries, no decorators, no explicit imports needed — just create a module
or package in `analyses/` that subclasses `Analysis`.

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

1. Run `polyzymd new-analysis <name>` to write the module and its tests
2. Declare `name`, `Settings`, `references` and `compute()` on a plain class
3. Return one `Observable` per reported quantity, each with a `kind` and a
   `unit`. The kind decides the reduction, the test and the figure
4. End the module with `NameAnalysis = contract_analysis(Name)`
5. **Test**: `PYTHONPATH=$PWD/src pixi run -e test pytest tests/analyses/plugins/test_<name>.py -q`
6. The CLI discovers it automatically through `polyzymd compare run <name>`

See `polyzymd.analyses.contract` for the contract,
`analysis-module.md` for detailed patterns, and
`docs/source/contributor_guide/extending_analyses.md` for the contributor
tutorial.

### Adding comparison statistics or formatters

Comparison statistics now live in the analyses framework itself. Use
`analyses/stats.py` for default scalar comparisons and
`analyses/shared/inferential_statistics.py` for reusable inferential helpers.
New plugins should keep formatting inline in the plugin's `format()` method.
Established plugins extract plotting into `_plotters.py` modules within each
plugin package; new plugins can start with plotting inline in `plot()` and
extract later as complexity grows.
