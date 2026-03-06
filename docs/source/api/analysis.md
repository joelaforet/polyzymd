# Analysis Module

Post-simulation analysis tools for enzyme-polymer MD trajectories.
Each analysis type has a **calculator** (computes the metric from a trajectory)
and a **result** class (serializable Pydantic model).

For tutorial-level documentation with worked examples, see the
{doc}`../tutorials/analysis_complete_workflow`.

## Core Infrastructure

Framework classes shared by all analyzers: the registry, base result model,
autocorrelation-aware statistics, and the MetricType system for handling
correlated MD data.

```{eval-rst}
.. automodule:: polyzymd.analysis.core.registry
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.core.metric_type
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.core.autocorrelation
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.core.statistics
   :members:
   :undoc-members:
   :show-inheritance:
```

## Results Base

```{eval-rst}
.. automodule:: polyzymd.analysis.results.base
   :members:
   :undoc-members:
   :show-inheritance:
```

## RMSF

Root Mean Square Fluctuation — per-residue structural deviation from a
catalytically competent reference.

See the {doc}`../tutorials/analysis_rmsf_quickstart` tutorial.

```{eval-rst}
.. automodule:: polyzymd.analysis.rmsf.calculator
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.results.rmsf
   :members:
   :undoc-members:
   :show-inheritance:
```

## Distances

Frame-by-frame inter-atomic or inter-group distances with KDE-based
distribution analysis.

See the {doc}`../tutorials/analysis_distances_quickstart` tutorial.

```{eval-rst}
.. automodule:: polyzymd.analysis.distances.calculator
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.results.distances
   :members:
   :undoc-members:
   :show-inheritance:
```

## Catalytic Triad

Simultaneous contact fraction for catalytic residue groups — measures how often
all pairwise distances are below threshold at the same time.

See the {doc}`../tutorials/analysis_triad_quickstart` tutorial.

```{eval-rst}
.. automodule:: polyzymd.analysis.triad.analyzer
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.results.triad
   :members:
   :undoc-members:
   :show-inheritance:
```

## Contacts

Polymer-protein contact detection with extensible criteria, per-residue
contact fractions, residence times, and interaction matrices.

See the {doc}`../tutorials/analysis_contacts_quickstart` tutorial.

```{eval-rst}
.. automodule:: polyzymd.analysis.contacts.calculator_parallel
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.contacts.calculator
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.contacts.results
   :members:
   :undoc-members:
   :show-inheritance:
```

```{eval-rst}
.. automodule:: polyzymd.analysis.contacts.aggregator
   :members:
   :undoc-members:
   :show-inheritance:
```

### Binding Preference

Enrichment ratios answering "does polymer type X preferentially bind amino acid
class Y?" — the foundation for the binding free energy and polymer affinity
comparisons.

See the {doc}`../tutorials/analysis_binding_preference` tutorial.

```{eval-rst}
.. automodule:: polyzymd.analysis.contacts.binding_preference
   :members:
   :undoc-members:
   :show-inheritance:
```

### Surface Exposure

Static surface exposure classification used to define the "expected" contact
share in enrichment calculations.

```{eval-rst}
.. automodule:: polyzymd.analysis.contacts.surface_exposure
   :members:
   :undoc-members:
   :show-inheritance:
```

### Contact Criteria

Strategy classes that define what constitutes a "contact" between a polymer
atom and a protein residue.

```{eval-rst}
.. automodule:: polyzymd.analysis.contacts.criteria
   :members:
   :undoc-members:
   :show-inheritance:
```

## SASA

Per-frame, per-residue Solvent Accessible Surface Area via the Shrake-Rupley
algorithm. Provides the exposure time series consumed by exposure dynamics.

```{eval-rst}
.. automodule:: polyzymd.analysis.sasa.trajectory
   :members:
   :undoc-members:
   :show-inheritance:
```

## Exposure Dynamics

Combines SASA time series with contact data to detect chaperone-like
polymer-protein interactions — residues that become transiently exposed and
are re-buried with the help of polymer contacts.

See the {doc}`../tutorials/analysis_exposure_dynamics` tutorial.

### Residue Classification

```{eval-rst}
.. automodule:: polyzymd.analysis.exposure.classification
   :members:
   :undoc-members:
   :show-inheritance:
```

### Dynamics Analysis

```{eval-rst}
.. automodule:: polyzymd.analysis.exposure.dynamics
   :members:
   :undoc-members:
   :show-inheritance:
```

### Chaperone Event Detection

```{eval-rst}
.. automodule:: polyzymd.analysis.exposure.chaperone
   :members:
   :undoc-members:
   :show-inheritance:
```

### Chaperone Enrichment

```{eval-rst}
.. automodule:: polyzymd.analysis.exposure.enrichment
   :members:
   :undoc-members:
   :show-inheritance:
```

## Common Utilities

Amino acid classification schemes and molecular selectors shared across
analysis types.

```{eval-rst}
.. automodule:: polyzymd.analysis.common.aa_classification
   :members:
   :undoc-members:
   :show-inheritance:
```

## Configuration

```{eval-rst}
.. automodule:: polyzymd.analysis.config
   :members:
   :undoc-members:
   :show-inheritance:
```
