# Documentation Rules

## Stack

- **Engine:** Sphinx 8.x
- **Parser:** MyST-Parser (Markdown support)
- **Theme:** sphinx_rtd_theme
- **Extensions:** autodoc, napoleon, intersphinx, myst-parser, copybutton
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
Incremental builds cause stale sidebar navigation (known issue).

## File Conventions

- New documentation files use MyST Markdown (`.md`)
- Existing `.rst` files (e.g., `index.rst`, `api/` stubs) can remain as-is
- Tutorials go in `docs/source/tutorials/`
- API reference goes in `docs/source/api/`
- Every new page must be added to a `toctree` directive

## Directory Structure

```
docs/source/
├── conf.py              # Sphinx configuration
├── index.rst            # Root toctree
├── api/                 # Auto-generated API docs
│   ├── index.rst
│   ├── analysis.md      # analysis/ package API
│   ├── analyses.md      # analyses/ plugin system API
│   └── *.rst            # Per-module API stubs
├── tutorials/           # User-facing guides
│   ├── architecture.md
│   ├── contributing.md
│   ├── extending_analyses.md    # ★ Primary guide for new contributors
│   ├── extending_comparators.md # Legacy (deprecated — see extending_analyses)
│   ├── extending_plotters.md    # Redirect — see extending_analyses
│   └── ...
├── contributor_guide/   # Contributor landing page
│   ├── index.md
│   └── setup.md
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
| `ReplicateContext` | `analyses/base.py` | Context for compute_replicate |
| `AggregateContext` | `analyses/base.py` | Context for aggregate |
| `ComparisonContext` | `analyses/base.py` | Context for compare |
| `PlotContext` | `analyses/base.py` | Context for plot |
| `MetricValue` | `analyses/base.py` | Scalar metric descriptor |
| `ComparisonResult` | `analyses/base.py` | Universal comparison result |
| `get_analysis()` | `analyses/discovery.py` | Plugin lookup |
| `list_analyses()` | `analyses/discovery.py` | Plugin enumeration |
| `run_comparison()` | `analyses/orchestrator.py` | Run one analysis |
| `run_all_comparisons()` | `analyses/orchestrator.py` | Run all analyses |

## ReadTheDocs

The project is configured for ReadTheDocs via `.readthedocs.yaml` at the repo
root. RTD builds use the pixi environment specification.
