# Analysis Shared Utilities

API reference for `polyzymd.analyses.shared` and its contributor-facing
submodules. Framework internals, CLI logging, and plugin-private artifact
helpers are documented with their owning packages.

## `polyzymd.analyses.shared`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.loader`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.loader
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.window`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.window
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.alignment`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.alignment
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.autocorrelation`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.autocorrelation
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.statistics`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.statistics
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.plotting`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.plotting
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.selections`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.selections
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Compatibility facades

`polyzymd.analyses.shared.config_hash` and `polyzymd.analyses.shared.sasa`
remain importable as temporary compatibility facades. New code should use the
owning modules instead:

- `polyzymd.analyses._framework.cache_identity` for framework cache identity
  helpers.
- `polyzymd.analyses.sasa._artifacts` for SASA plugin artifact helpers.

## `polyzymd.analyses.shared.groupings`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.groupings
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## `polyzymd.analyses.shared.selectors`

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.selectors
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
