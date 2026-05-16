# Documentation Rules

## Stack

- **Engine:** Sphinx 8.x
- **Parser:** MyST-Parser (Markdown support)
- **Theme:** sphinx_rtd_theme
- **Extensions:** autodoc, autodoc-pydantic, napoleon, intersphinx, myst-parser,
  copybutton, sphinx-design, autosummary, viewcode, mathjax
- **Pydantic rendering:** `sphinxcontrib.autodoc_pydantic` handles Pydantic v2
  `BaseModel` field documentation
- **Source:** `docs/source/`
- **Build output:** `docs/build/html/`
- **Config:** `docs/source/conf.py`

## Building Docs

```bash
# Full clean build (always use this after structural changes)
pixi run -e build make -C docs clean html

# Quick rebuild (only for content-only changes, NOT structural)
pixi run -e build make -C docs html
```

**Important:** After adding new pages to toctree, ALWAYS use `make clean html`.
Incremental builds cause stale sidebar navigation.

## Zero-Warning Build Gate

The doc build **must produce zero warnings**. After every change, verify:

```bash
pixi run -e build make -C docs clean html 2>&1 | grep -E "(WARNING|ERROR|build succeeded)"
```

Output should show only `build succeeded.` with no warning count. If you see
`build succeeded, N warnings.` all N warnings must be fixed.

## Warning Prevention Rules

### Duplicate Object Description Warnings

The most common source of warnings. Sphinx autodoc generates duplicate
`attribute` directives when a module contains Pydantic `BaseModel` subclasses
or `@dataclass` classes.

**Root cause:** `autodoc_typehints = "description"` combined with
`special-members: __init__` generates attribute entries from `__init__`
parameters, which collide with class-body field documentation (or
autodoc-pydantic's `pydantic_field` directives for Pydantic models).

**Fix:** Add `:no-index:` to the `automodule` directive on the API page:

    ```{eval-rst}
    .. automodule:: polyzymd.some.module
       :members:
       :undoc-members:
       :show-inheritance:
       :no-index:
    ```

**Modules currently using `:no-index:`:**

| API Page | Module | Reason |
|----------|--------|--------|
| `api/config.md` | `polyzymd.config.schema` | 24 Pydantic BaseModel classes |
| `api/workflow.md` | `polyzymd.workflow.slurm` | SlurmConfig, JobContext dataclasses |
| `api/workflow.md` | `polyzymd.workflow.daisy_chain` | DaisyChainConfig, SubmissionResult dataclasses |
| `api/workflow.md` | `polyzymd.workflow.analysis_slurm` | Dataclasses for analysis SLURM jobs |
| `api/core.md` | `polyzymd.core.parameters` | SimulationParameters and related dataclasses |
| `api/core.md` | `polyzymd.core.restraints` | RestraintDefinition, AtomSelection dataclasses |
| `api/builders.md` | `polyzymd.builders.solvent` | SolventComposition, SolvationCounts dataclasses |

**Rule:** When adding a new `automodule` directive, check whether the module
contains Pydantic `BaseModel` subclasses or `@dataclass` classes. If yes, add
`:no-index:`.

### Non-Consecutive Header Levels

MyST Markdown enforces consecutive heading levels. Going from `##` to `####`
(skipping `###`) produces a warning. Always use consecutive levels.

### Docstring reST Formatting (Prevents "Unexpected indentation" Errors)

Python docstrings are rendered as reStructuredText by autodoc. Follow these
rules:

1. **Blank line before lists:** Always insert a blank line between a paragraph
   and a bulleted or enumerated list.
2. **Consistent sub-bullet indentation:** Sub-bullets align to the parent
   item's text start position.
3. **Continuation lines:** Wrapped lines under numbered items align to the
   first character after `N. `.
4. **Pydantic model docstrings:** Do not manually document fields in
   `Attributes:` sections — autodoc-pydantic handles field documentation
   automatically.

### Autodoc Import Failures

If a new dependency is used at module level and isn't in `autodoc_mock_imports`
in `conf.py`, the doc build will fail. Add it to the mock list.

## File Conventions

- New documentation files use MyST Markdown (`.md`)
- Existing `.rst` files (e.g., `index.rst`, `api/` stubs) can remain as-is
- Every new page must be added to a `toctree` directive

## Directory Structure

```
docs/source/
├── conf.py              # Sphinx configuration
├── index.rst            # Root toctree
├── get_started/         # Installation, quickstart
├── tutorials/           # Learning-oriented guided walkthroughs
├── how_to/              # Task-oriented practical guides
├── reference/           # Lookup-oriented factual docs
│   ├── data_requirements.md    # Data & directory layout reference
│   ├── comparison_yaml.md      # comparison.yaml schema reference
│   └── cli_reference.md        # CLI command reference
├── explanation/         # Conceptual "why" discussions
├── contributor_guide/   # Contributor landing page and setup
├── api/                 # Auto-generated API docs
│   ├── index.rst
│   ├── config.md        # config/ (has :no-index: for Pydantic models)
│   ├── workflow.md      # workflow/ (has :no-index: for dataclasses)
│   ├── core.md          # core/ (has :no-index: for dataclasses)
│   ├── builders.md      # builders/ (solvent has :no-index:)
│   ├── analyses.md      # analyses/ plugin system
│   └── *.rst            # Per-module stubs
└── _static/             # Static assets (CSS, images)
```

## Docstring Style

NumPy-style for all new/updated docstrings. Napoleon extension handles
rendering in Sphinx. See `code-style.md` for the full docstring template.

## Cross-References

Use MyST roles for cross-referencing:
- `{doc}path/to/page` — link to another page
- `{ref}label` — link to a labeled section
- `{func}polyzymd.analyses.discovery.get_analysis` — link to API docs
- `{class}polyzymd.analyses.base.Analysis` — link to class docs

## Key API Classes for Documentation

| Class | Location | Role |
|-------|----------|------|
| `Analysis` | `analyses/base.py` | Plugin base class |
| `ReplicateContext` | `analyses/base.py` | Context for per-replicate execution |
| `AggregateContext` | `analyses/base.py` | Context for aggregate |
| `ComparisonContext` | `analyses/base.py` | Context for compare |
| `PlotContext` | `analyses/base.py` | Context for plot |
| `MetricValue` | `analyses/base.py` | Scalar metric descriptor |
| `ComparisonResult` | `analyses/base.py` | Universal comparison result |
| `get_analysis()` | `analyses/discovery.py` | Plugin lookup |
| `list_analyses()` | `analyses/discovery.py` | Plugin enumeration |
| `run_comparison()` | `analyses/orchestrator.py` | Run one analysis |
| `run_all_comparisons()` | `analyses/orchestrator.py` | Run all analyses |

`polyzymd.analyses.base` is a public facade. It re-exports context and result
models from private framework modules so contributor documentation can keep a
single stable import path. Do not document `_analysis_*`, `_contexts.py`, or
`_comparison_models.py` as user-facing API pages unless explicitly writing
internal developer reference material.

## ReadTheDocs

The project is configured for ReadTheDocs via `.readthedocs.yaml` at the repo
root. RTD builds use the pixi environment specification.
