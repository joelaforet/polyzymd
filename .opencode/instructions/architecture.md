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
| **Plugins** (public) | `rmsf/`, `contacts/`, `distances/`, etc. | One class per analysis type — the extension point |
| **Private modules** | `_analysis_*`, `_contexts.py`, `_comparison_models.py`, `<name>/_*.py`, etc. | Internal framework and plugin implementation details; not contributor import targets |
| **Shared utilities** | `shared/loader.py`, `shared/alignment.py`, etc. | `TrajectoryLoader`, alignment, statistics |
| **Framework** | `base.py`, `discovery.py`, `orchestrator.py`, `stats.py`, `mda/` | Stable public facade, auto-discovery, artifact lifecycle |

New analysis types may be simple single-file modules or packages under
`analyses/`. Compute-stage plugins use `build_mda_jobs()` to create
`MDAAnalysisJob` objects around `AnalysisBase`-compatible work and, when
needed, `build_mda_collector()` to map completed jobs to `ReplicateArtifact`.
PolyzyMD owns `ArtifactStore`, `ConditionArtifact`, `ComparisonArtifact`,
ensemble orchestration, statistics, and plotting. Advanced packages should keep
MDAnalysis helpers in a dedicated module such as `_mda.py`.

`polyzymd.analyses.base` is the stable public API facade for contributors. It
re-exports `Analysis`, lifecycle context objects, metric descriptors, and
comparison models while delegating implementation to private modules such as
`_analysis_compare.py`, `_analysis_io.py`,
`_analysis_contract.py`, `_contexts.py`, and `_comparison_models.py`. Do not
import these private modules from contributor plugins; import public symbols
from `polyzymd.analyses.base`.

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

1. Run `polyzymd new-analysis <name>` to scaffold the plugin, OR create a module/package under `src/polyzymd/analyses/` manually
2. Subclass `Analysis` with `name` and `Settings`
3. Choose the lifecycle mode:
   - MDAnalysis-native plugin: implement `build_mda_jobs()` and, when needed, `build_mda_collector()` when `has_compute_stage=True`; MDAnalysis owns per-trajectory iteration while PolyzyMD owns artifacts, ensemble aggregation, statistics, and plotting
   - Compare-only plugin: set `has_compute_stage=False`
4. Implement `aggregate()` only when `has_aggregate_stage=True`
5. For default comparison: implement `extract_metrics()`
6. Implement `plot()` using `_build_plot_data()` and shared plotting helpers
7. Optionally implement `format()`
8. **Test**: `pixi run -e build pytest tests/analyses/plugins/test_<name>.py -v`
9. The CLI automatically discovers it via `polyzymd compare run <name>`

See the `polyzymd.analyses.base` facade for the public contract,
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
