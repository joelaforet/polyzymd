# Experimental analyses

PolyzyMD v1.3 ships the stable analyses documented in the analysis reference
pages. The following experimental analyses were removed from the shipped v1.3
code and are not available:

- contacts binding preference
- exposure dynamics
- binding free energy
- polymer affinity
- polymer bridging

These archived features are retained only for historical review. They are not
maintained as part of the v1.3 analysis workflow, and this page intentionally
does not reproduce the old how-to content.

This is an intentional breaking change for v1.3. The archived analyses were
experimental, had no known external users, and did not yet have stable result or
cache contracts. PolyzyMD therefore does not migrate their old configuration
blocks, cached intermediates, or result files forward into the active analysis
workflow.

## Archived names and migration path

If you need a related stable analysis, recompute from
the original trajectories with the analysis listed below; do not reuse old
experimental caches as input to v1.3 comparisons.

| Archived request names | Active replacement path |
|---|---|
| `binding_preference`, `contacts_binding_preference`, `contact_binding_preference`, `bp` | No one-to-one replacement. Run `polyzymd analyze contacts`, whose `<type>_contact_fraction` results give the contact fraction with each monomer type, and perform project-specific post-processing if you need binding-preference summaries. |
| `exposure`, `exposure_dynamics`, `surface_exposure` | No direct replacement. Run `polyzymd analyze sasa` for solvent exposure observables when SASA answers the scientific question. |
| `binding_free_energy`, `bfe` | No active PolyzyMD replacement. Use a dedicated free-energy workflow outside v1.3, or inspect the archived implementation for historical context only. |
| `polymer_affinity`, `pa` | No one-to-one replacement. Run `polyzymd analyze contacts` and/or `polyzymd analyze distances` when those observables are sufficient. |
| `polymer_bridging`, `bridging` | No one-to-one replacement. Run `polyzymd analyze contacts` and/or `polyzymd analyze distances` and implement bridging-specific logic externally if required. |

If you previously ran one of the archived analyses, keep those outputs as
historical artifacts. For v1.3 analyses, clear any experimental cache
directories of those analyses and rerun the active
replacement analysis from trajectories so the result metadata, cache identity,
and comparison statistics all come from the current code.

## Current workflow

Every shipped analysis, hydrogen bonds included, runs through
`polyzymd analyze <analysis> -c <config.yaml> ...`, which draws its figures;
see {doc}`analysis_functions`. `polyzymd analyze` refuses the archived names as
it refuses any unknown analysis name, listing the names it runs.
