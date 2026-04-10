# Explanation

Explanation pages help you understand why PolyzyMD works the way it does, how
to interpret outputs, and what assumptions or tradeoffs matter.

## Concepts and Design

- [Architecture](architecture.md)
- [Residue Assignment and Chain Conventions](residue_assignment.md)
- [Colored Logging](colored_logging.md)

## Analysis Interpretation and Best Practices

- [Statistical Best Practices for Analysis](analysis_statistics_best_practices.md)
- [RMSD Best Practices](analysis_rmsd_best_practices.md)
- [Establishing Convergence in MD Simulations](convergence_detection.md) — Conceptual guide to convergence detection: what it means, how PolyzyMD's sliding-window heuristic works, and when to tune parameters.
- [Rg Best Practices](analysis_rg_best_practices.md)
- [RMSF Best Practices](analysis_rmsf_best_practices.md)
- [RMSF Reference Selection](analysis_reference_selection.md)
- [RMSF Verification](analysis_rmsf_verification.md)
- [Catalytic Triad Best Practices](analysis_triad_best_practices.md)

## Experimental Analysis Concepts

These workflows remain available in PolyzyMD, but their interpretation is still
evolving and they should be treated as experimental:

- [Binding Preference](../how_to/analysis_binding_preference.md)
- [Binding Free Energy](../how_to/analysis_binding_free_energy.md)
- [Polymer Affinity](../how_to/analysis_polymer_affinity.md)
- [Polymer Bridging Interpretation](polymer_bridging_interpretation.md)
- [Exposure Dynamics](../how_to/analysis_exposure_dynamics.md)

<!-- IMAGE OPPORTUNITY: Add a chain-convention figure (A/B/C/D+), plus a simple
stable-vs-experimental analysis map to orient users before they dive into the
longer conceptual pages. -->

```{toctree}
:hidden:
:maxdepth: 1

Architecture <architecture>
Residue Assignment and Chain Conventions <residue_assignment>
Colored Logging <colored_logging>
Statistical Best Practices for Analysis <analysis_statistics_best_practices>
RMSD Best Practices <analysis_rmsd_best_practices>
Convergence Detection <convergence_detection>
Rg Best Practices <analysis_rg_best_practices>
RMSF Best Practices <analysis_rmsf_best_practices>
RMSF Reference Selection <analysis_reference_selection>
RMSF Verification <analysis_rmsf_verification>
Catalytic Triad Best Practices <analysis_triad_best_practices>
Polymer Bridging Interpretation <polymer_bridging_interpretation>
```
