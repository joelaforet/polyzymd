"""Contact analysis module for polymer-protein interactions.

This module provides comprehensive contact analysis between:
- Polymer residues and protein residues
- Protein residues and solvent/cosolvent
- Substrate/ligand and protein or polymer

Key features:
- Strategy pattern for contact criteria (easily extensible)
- Compressed event storage (start_frame, duration)
- Frame-by-frame computation for proper autocorrelation analysis
- Multi-replicate aggregation with proper statistics
- **Parallel analysis** with optimized neighbor searching
- **Statistical analysis** following LiveCoMS best practices (Grossfield et al. 2018)
- **Autocorrelation-corrected** uncertainty quantification

Topology and Chain Conventions
------------------------------
PolyzyMD uses a consistent chain naming convention:

==========  =====================================
Chain ID    Contents
==========  =====================================
Chain A     Protein/Enzyme
Chain B     Substrate/Ligand
Chain C     Polymers
Chain D+    Solvent (water, ions, co-solvents)
==========  =====================================

**IMPORTANT**: Always use ``solvated_system.pdb`` as the topology file for
contact analysis, NOT ``production_N_topology.pdb``. The solvated_system.pdb
contains correct chain assignments from system setup, while production
topologies may have modified or lost chain information.

Example topology loading:

>>> import MDAnalysis as mda
>>> topology = project_dir / "solvated_system.pdb"  # Correct
>>> # topology = production_dir / "production_0_topology.pdb"  # WRONG!
>>> u = mda.Universe(str(topology), trajectory_files)

Default selections based on chain convention:

>>> polymer_sel = "segid C"          # All polymer atoms (Chain C)
>>> protein_sel = "protein"          # All protein atoms (usually Chain A)
>>> aromatic_sel = "protein and (resname TRP PHE TYR HIS)"  # Aromatic residues

Statistical Analysis
--------------------
Contact analysis includes proper uncertainty quantification following
LiveCoMS best practices (Grossfield et al. 2018):

- **Statistical inefficiency (g)**: Computed via autocorrelation function
  integration with finite-size correction (Chodera et al. 2007)
- **Effective sample size**: N_eff = N/g accounts for correlation
- **Reliability warning**: Warns if N_eff < 10 (results may be unreliable)
- **Residence time statistics**: Mean, SEM, and distribution of contact durations

Examples
--------
>>> from polyzymd.analysis.contacts import ContactAnalyzer, AnyAtomWithinCutoff
>>> from polyzymd.analysis.common.selectors import ProteinResidues, PolymerChains
>>>
>>> # Create analyzer with default criteria
>>> analyzer = ContactAnalyzer(
...     target_selector=ProteinResidues(),
...     query_selector=PolymerChains(),
...     criteria=AnyAtomWithinCutoff(cutoff=4.0),
... )
>>>
>>> # Run analysis
>>> result = analyzer.run(universe)
>>> print(f"Total contacts: {result.n_contact_events}")

For large systems, use the parallel analyzer:

>>> from polyzymd.analysis.contacts import ParallelContactAnalyzer
>>>
>>> analyzer = ParallelContactAnalyzer(
...     target_selector=ProteinResidues(),
...     query_selector=PolymerChains(),
...     cutoff=4.0,
... )
>>>
>>> # Run with 4 parallel workers (~4x speedup)
>>> result = analyzer.run(universe, n_jobs=4)

For multi-replicate aggregation with statistics:

>>> from polyzymd.analysis.contacts.aggregator import aggregate_contact_results
>>>
>>> # Aggregate results from multiple replicates
>>> agg = aggregate_contact_results(results, replicates=[1, 2, 3])
>>> print(f"Contact fraction: {agg.contact_fraction_mean:.1%} ± {agg.contact_fraction_sem:.1%}")
>>> print(f"Statistical inefficiency (g): {agg.statistical_inefficiency:.2f}")
>>> print(f"Effective samples: {agg.n_effective:.1f}")

References
----------
.. [1] Grossfield, A. et al. (2018). Best practices for quantification of
       uncertainty and sampling quality in molecular simulations. *Living
       Journal of Computational Molecular Science*, 1(1), 5067.
       https://doi.org/10.33011/livecoms.1.1.5067

.. [2] Chodera, J. D. et al. (2007). Use of the weighted histogram analysis
       method for the analysis of simulated and parallel tempering simulations.
       *Journal of Chemical Theory and Computation*, 3(1), 26-41.
       https://doi.org/10.1021/ct0502864
"""

from polyzymd.analysis.contacts.criteria import (
    ContactCriteria,
    AnyAtomWithinCutoff,
    AnyAtomToCOM,
    COMToCOM,
    MinimumDistance,
)
from polyzymd.analysis.contacts.results import (
    ContactEvent,
    ResidueContactData,
    ContactResult,
)
from polyzymd.analysis.contacts.calculator import ContactAnalyzer
from polyzymd.analysis.contacts.calculator_parallel import ParallelContactAnalyzer
from polyzymd.analysis.contacts.surface_exposure import (
    ResidueExposure,
    SurfaceExposureResult,
    SurfaceExposureFilter,
)
from polyzymd.analysis.contacts.binding_preference import (
    BindingPreferenceEntry,
    BindingPreferenceResult,
    AggregatedBindingPreferenceEntry,
    AggregatedBindingPreferenceResult,
    compute_binding_preference,
    aggregate_binding_preference,
    resolve_protein_group_selections,
    resolve_polymer_type_selections,
    compute_binding_preference_from_config,
)

__all__ = [
    # Criteria (Strategy pattern)
    "ContactCriteria",
    "AnyAtomWithinCutoff",
    "AnyAtomToCOM",
    "COMToCOM",
    "MinimumDistance",
    # Results
    "ContactEvent",
    "ResidueContactData",
    "ContactResult",
    # Calculators
    "ContactAnalyzer",
    "ParallelContactAnalyzer",
    # Surface Exposure
    "ResidueExposure",
    "SurfaceExposureResult",
    "SurfaceExposureFilter",
    # Binding Preference
    "BindingPreferenceEntry",
    "BindingPreferenceResult",
    "AggregatedBindingPreferenceEntry",
    "AggregatedBindingPreferenceResult",
    "compute_binding_preference",
    "aggregate_binding_preference",
    "resolve_protein_group_selections",
    "resolve_polymer_type_selections",
    "compute_binding_preference_from_config",
]
