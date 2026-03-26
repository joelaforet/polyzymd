# Code Style Rules

## Formatting

- **Black** for formatting, line-length=100
- **Ruff** for linting (see `pyproject.toml` `[tool.ruff]` section for enabled rules)
- Run `ruff check src/` and `black src/ --check` before committing
- Auto-fix with `ruff check src/ --fix` and `black src/`

## Import Conventions

Order: stdlib, third-party, local (separated by blank lines).

```python
# stdlib
import logging
from pathlib import Path
from typing import Any

# third-party
import numpy as np
import pandas as pd

# local
from polyzymd.config.schema import SimulationConfig
from polyzymd.core.base import BaseAnalyzer
```

**Lazy imports** are mandatory for heavy dependencies. These packages must NOT
appear at module level:

- `openmm` / `openmm.app` / `openmm.unit`
- `openff.toolkit` / `openff.interchange`
- `MDAnalysis` (`mda`)
- `parmed`
- `pdbfixer`

```python
# CORRECT - lazy import inside function
def run_simulation(config):
    from openmm import app, unit
    from openmm.app import Simulation
    ...

# WRONG - module-level import of heavy dep
import openmm  # This breaks environments without OpenMM installed
```

## Type Annotations

- Use Python 3.10+ union syntax: `X | None` instead of `Optional[X]`
- Use `from __future__ import annotations` if needed for forward references
- Older modules still use `Optional[X]` — don't refactor unless touching that code

## Naming

- Classes: `PascalCase` (e.g., `ContactAnalyzer`, `SystemBuilder`)
- Functions/methods: `snake_case` (e.g., `calculate_contacts`, `from_config`)
- Constants: `UPPER_SNAKE_CASE` (e.g., `DEFAULT_CUTOFF`, `CHAIN_PROTEIN`)
- Private: prefix with underscore (e.g., `_validate_input`, `_cache`)

## Docstrings

**NumPy style** for all new code:

```python
def example(param_a: str, param_b: int = 10) -> bool:
    """Short summary on one line.

    Extended description if needed, explaining behavior,
    edge cases, or algorithm details.

    Parameters
    ----------
    param_a : str
        Description of param_a.
    param_b : int, optional
        Description of param_b, by default 10.

    Returns
    -------
    bool
        Description of return value.
    """
```

Google-style docstrings exist in older modules (`builders/`, `simulation/`,
`workflow/`). Don't convert them unless you're significantly modifying the
function.
