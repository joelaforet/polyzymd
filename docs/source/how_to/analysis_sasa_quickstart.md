# SASA Analysis: Quick Start

Use the `sasa` plugin to measure solvent-accessible surface area for whole
proteins, active sites, or polymer-shielded regions.

```{note}
For a guided learning path, see {doc}`../tutorials/sasa_analysis`. For settings,
artifact paths, and output fields, see {doc}`../reference/analysis_sasa_reference`.
```

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## TL;DR

```bash
# Run only SASA for conditions in comparison.yaml
polyzymd compare run sasa -f comparison.yaml --eq-time 10ns

# Run all enabled analyses, including SASA
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recomputation when settings or selections changed
polyzymd compare run sasa -f comparison.yaml --eq-time 10ns --recompute
```

## Prerequisites

Before running SASA, confirm you have:

1. completed production trajectories,
2. a `comparison.yaml` with one or more conditions,
3. valid simulation `config.yaml` paths for each condition, and
4. selections that match your topology.

PolyzyMD examples use the chain convention A = protein, B = substrate,
C = polymer, and D+ = solvent/ions/other.

## Configure a minimal whole-protein run

Use this when you only need total protein SASA.

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_total"
        target_selection: "protein"
```

When `context_selection` is omitted, it defaults to the same value as
`target_selection`. This reports the protein's self-SASA.

Run it:

```bash
polyzymd compare run sasa -f comparison.yaml --eq-time 10ns
```

## Configure a two-run shielding comparison

Use this when you want a practical polymer-shielding signal.

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"
      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainid C"
```

Interpretation:

- `protein_isolated` is the baseline protein surface.
- `protein_with_polymer` allows polymer atoms to block protein surface points.
- A lower `protein_with_polymer` value in polymer conditions indicates
  shielding.

## Focus on active-site exposure

Use residue selections when the biological question is whether polymer blocks a
catalytic site or binding pocket.

```yaml
plugins:
  sasa:
    runs:
      - label: "active_site_isolated"
        target_selection: "protein and (resid 77 or resid 156 or resid 262)"
        context_selection: "protein"
      - label: "active_site_with_polymer"
        target_selection: "protein and (resid 77 or resid 156 or resid 262)"
        context_selection: "protein or chainid C"
```

Adjust residue IDs to match your enzyme. If the polymer-aware active-site run
has lower SASA, polymer may be reducing access to that site.

## Compare monomer-specific shielding

Use monomer residue names when your polymer contains distinct monomer types.

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"
      - label: "protein_with_sbma"
        target_selection: "protein"
        context_selection: "protein or resname SBMA"
      - label: "protein_with_egma"
        target_selection: "protein"
        context_selection: "protein or resname EGMA"
```

Check topology residue names before relying on a monomer-specific selection:

```bash
python - <<'PY'
import MDAnalysis as mda

u = mda.Universe("solvated_system.pdb")
print(sorted(set(u.select_atoms("chainid C").residues.resnames)))
PY
```

## Keep memory bounded on long trajectories

SASA is CPU-intensive. `chunk_size` sets how many frames go to MDTraj at once,
which bounds memory. It is not free of consequence: MDTraj returns slightly
different areas for the same frame depending on the size of the array it
arrives in, worth about 0.1 percent of the total, so use one value for every
condition you intend to compare.

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainid C"
    chunk_size: 50
    n_sphere_points: 960
```

Practical guidance:

- Lower `chunk_size` if memory is tight, and then use the same value for every
  condition of the comparison.
- Lower `n_sphere_points` for an exploratory scan, then restore 960 for the
  final numbers.
- Narrow the window with `--eq-time` rather than a per-run stride. A per-run
  `stride` is deprecated and ignored, so every observable shares one frame
  window and the conditions stay comparable.

## Run on SLURM instead of locally

For large systems or many replicates, submit analysis jobs to SLURM:

```bash
polyzymd compare submit sasa -f comparison.yaml --dry-run
```

Inspect the generated jobs, then submit without `--dry-run` when the resource
requests look right. SASA has a high execution-cost hint, so use the full HPC
guide for scheduler options, monitoring, and troubleshooting:
{doc}`hpc_execution`.

## Know why there are no plot files yet

Figures are drawn from the observable kind in the contract runner, which is
still being written. A `sasa` run today writes artifacts and the text report but
no figures. A `plot_settings.sasa` block in an existing comparison file still
loads and warns that it does nothing.

## Quick output checks

After a run, confirm the canonical outputs exist:

```bash
ls analysis/<condition>/sasa/run_1/
ls analysis/<condition>/sasa/aggregated/
ls comparison/sasa/
```

Inspect the condition aggregates from the comparison result:

```bash
python - <<'PY'
import json
from pathlib import Path

result = json.loads(Path("comparison/sasa/result.json").read_text())
for label, observables in result["payload"]["conditions"].items():
    print(label)
    for observable in observables:
        if observable["kind"] == "mean_of_timeseries":
            print(
                f"  {observable['name']}: {observable['mean']:.1f}"
                f" +/- {observable['sem']:.1f} {observable['unit']}"
                f" (n={observable['n_replicates']} replicates)"
            )
PY
```

Check the pairwise tests:

```bash
python - <<'PY'
import json
from pathlib import Path

result = json.loads(Path("comparison/sasa/result.json").read_text())
for comparison in result["payload"]["comparisons"]:
    print(
        comparison["name"],
        comparison["control"], "vs", comparison["condition"],
        f"delta {comparison['delta']:+.1f} {comparison['unit']}",
        f"p_adj {comparison['p_adjusted']:.3g}",
    )
PY
```

## Common fixes

| Symptom | Fix |
|---------|-----|
| A run raises `SelectionError` for zero atoms | Test the `target_selection` and `context_selection` against the topology. |
| SASA is too slow | Shorten the window with `--eq-time`, lower `n_sphere_points` for exploration, or submit to SLURM. |
| Memory use is too high | Lower `chunk_size`. |
| A run raises `ReplicateError` naming Tien | The target holds a residue with no maximum accessible area; restrict `target_selection` to standard amino acids. |
| Monomer-specific run looks empty | Verify `resname` values and chain C membership in the topology. |
| Changed selections but results did not change | Re-run with `--recompute`. |

## Where to find details

- Guided shielding tutorial: {doc}`../tutorials/sasa_analysis`
- Settings and output reference: {doc}`../reference/analysis_sasa_reference`
- Comparison file setup: {doc}`analysis_compare_conditions`
- SLURM execution: {doc}`hpc_execution`
