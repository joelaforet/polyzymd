# How to Compare Simulation Conditions

Use this guide when you already have completed PolyzyMD simulations and want to
run the current `polyzymd compare` workflow.

You will:

- create a comparison workspace
- configure one or more analysis plugins under `plugins:`
- run `polyzymd compare run` or `polyzymd compare run-all`
- generate figures with `polyzymd compare plot-all`

```{important}
For the `v1.3.0` release, the stable comparison stack is RMSD, Rg, RMSF,
contacts, distances, catalytic triad, secondary structure, and SASA. Binding
preference, exposure dynamics, binding free energy, and polymer affinity remain
available, but PolyzyMD labels them as experimental.
```

```{note}
If you have not yet run a full analysis/comparison workflow, start with
[Tutorial: Analyze a Study from Finished Simulations](../tutorials/analysis_complete_workflow.md).
```

:::{admonition} Environment Setup
:class: tip

All commands below assume you have activated the PolyzyMD pixi environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

## Before You Start

Make sure each condition already has:

- a simulation `config.yaml`
- finished trajectories for the replicates you want to compare
- any shared inputs needed by the plugin you plan to run

The comparison pipeline can reuse cached analysis data when it exists, but it
can also compute missing per-condition results during `polyzymd compare run`.

## Step 1: Create a Comparison Workspace

```bash
polyzymd compare init -n polymer_stability_study
cd polymer_stability_study
```

This creates:

```text
polymer_stability_study/
├── comparison.yaml
├── comparison/
├── figures/
└── structures/
```

- `comparison.yaml` defines the conditions and enabled plugins
- `comparison/` stores cached comparison JSON, one subdirectory per analysis
- `figures/` stores generated plots
- `structures/` holds shared reference files such as an enzyme PDB for SASA

## Step 2: Define a Minimal `comparison.yaml`

Start with one stable analysis. RMSF is a good first comparison because it has
few extra inputs.

```yaml
name: "polymer_stability_study"
description: "Effect of polymer composition on enzyme flexibility"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_enzyme_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../SBMA_100_enzyme_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../EGMA_100_enzyme_DMSO/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
```

To enable more analyses, add more sections under `plugins:`:

```yaml
plugins:
  rmsf:
    selection: "protein and name CA"

  contacts:
    polymer_selection: "chainID C"
    protein_selection: "protein"
    cutoff: 4.5

  catalytic_triad:
    name: "Ser-His-Asp"
    threshold: 3.5
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"

  distances:
    pairs:
      - label: "Substrate-Ser77"
        selection_a: "resname SUB and name C1"
        selection_b: "protein and resid 77 and name OG"

  rmsd:
    runs:
      - label: "Protein Backbone"
        selection: "protein and name CA"
        alignment_selection: "protein and name CA"
        reference_mode: "centroid"
      - label: "Active Site"
        selection: "protein and (resid 77 or resid 133 or resid 156) and name CA"
        alignment_selection: "protein and name CA"
        reference_mode: "centroid"

  rg:
    runs:
      - label: "Whole Protein"
        selection: "protein"
      - label: "Protein Backbone"
        selection: "protein and name CA"
```

:::{admonition} Statistical settings for pairwise comparisons
:class: tip

Plugins that perform cross-condition statistical tests support per-plugin
settings in the `plugins:` block. For example, contacts supports `fdr_alpha`,
`min_effect_size`, and `top_residues`; binding free energy and polymer affinity
support `fdr_alpha`. See the
[Comparison Reference](../reference/analysis_comparison_reference.md#per-plugin-statistical-settings)
for the full settings table. For post-hoc method details (BH t-tests, Tukey
HSD, Cohen's d, and significance markers), see the
[Post-Hoc Testing Reference](../reference/posthoc_testing.md).
:::

## Step 3: Validate the Config

```bash
polyzymd compare validate
```

You should see a passing summary with the study name, condition count, and the
enabled plugin sections.

## Step 4: Run One Comparison

```bash
polyzymd compare run rmsf
```

This command:

- resolves `plugins.rmsf` from `comparison.yaml`
- computes or reloads per-condition RMSF data
- performs the cross-condition comparison
- writes the canonical cache file to `comparison/rmsf/result.json`
- prints a formatted summary to the terminal

:::{admonition} Running on an HPC cluster?
:class: tip

For expensive analyses (SASA, contacts, hydrogen bonds) or large studies with
many conditions and replicates, use `polyzymd compare submit` to dispatch
analysis as SLURM jobs instead of running interactively:

```bash
polyzymd compare submit sasa --partition <part> --mem 8G --time 02:00:00
polyzymd compare status sasa       # monitor progress
polyzymd compare finalize sasa     # (if needed) re-run compare + plot
```

Each replicate runs as an independent job, with automatic dependency wiring
for aggregation and finalization. See {doc}`hpc_execution` for the full
workflow, including dry-run previews and job arrays.
:::

You can save the formatted report separately with `-o`:

```bash
polyzymd compare run rmsf --format markdown -o reports/rmsf.md
```

## Step 5: Run All Enabled Comparisons

Once you have multiple plugin sections configured, run them together:

```bash
polyzymd compare run-all
```

Or run them and generate plots in one pass:

```bash
polyzymd compare run-all --plot
```

## Step 6: Generate Figures

For a plotting smoke test:

```bash
polyzymd compare plot-all --list-available
polyzymd compare plot-all
```

`--list-available` is useful because it shows which plot types are available
for the currently enabled plugins and which are experimental.

## Step 7: Check the Outputs

After a successful run, expect files like these:

```text
polymer_stability_study/
├── comparison.yaml
├── comparison/
│   ├── rmsf/
│   │   └── result.json
│   ├── contacts/
│   │   └── result.json
│   ├── distances/
│   │   └── result.json
│   └── catalytic_triad/
│       └── result.json
└── figures/
    ├── rmsf/
    │   ├── rmsf_comparison.png
    │   └── rmsf_profile.png
    └── ...
```

If your smoke test is `polyzymd compare plot-all`, success means:

- the command completes without error
- stable plots render normally
- experimental plots, if enabled, render with explicit experimental labeling

## Programmatic Use

If you need to run the comparison pipeline from Python, use the plugin
orchestrator directly:

```python
from pathlib import Path

from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.orchestrator import run_comparison
from polyzymd.config.comparison import ComparisonConfig

config = ComparisonConfig.from_yaml(Path("comparison.yaml"))
analysis = get_analysis("rmsf")()

pipeline_result = run_comparison(
    analysis,
    config,
    equilibration="10ns",
)

result = pipeline_result["comparison"]
print(result.ranking)
print(pipeline_result["comparison_path"])
```

## Adding More Stable Analyses

Common next additions to `comparison.yaml` are:

- `rmsd` for RMSD timeseries and structural stability comparison
- `rg` for Radius of Gyration and structural compactness comparison
- `contacts` for polymer coverage and contact fraction
- `distances` for custom atom-pair distances
- `catalytic_triad` for active-site geometry
- `secondary_structure` for helix/strand persistence and content

For end-to-end examples, see:

- [Run RMSD Analysis](analysis_rmsd_quickstart.md)
- [Run Rg Analysis](analysis_rg_quickstart.md)
- [Run RMSF Analysis](analysis_rmsf_quickstart.md)
- [Run Contacts Analysis](analysis_contacts_quickstart.md)
- [Run Distance Analysis](analysis_distances_quickstart.md)
- [Run Catalytic Triad Analysis](analysis_triad_quickstart.md)

## Experimental Workflows

Experimental workflows remain available, but they are not the default path for
the presentation release:

- [Experimental: Analyze Binding Preference](analysis_binding_preference.md)
- [Experimental: Analyze Binding Free Energy](analysis_binding_free_energy.md)
- [Experimental: Analyze Polymer Affinity](analysis_polymer_affinity.md)
- [Experimental: Analyze Polymer Bridging](analysis_polymer_bridging.md)
- [Experimental: Analyze Exposure Dynamics](analysis_exposure_dynamics.md)

## Troubleshooting

### `config` path not found

Paths in `comparison.yaml` are resolved relative to the location of
`comparison.yaml`, not your current shell directory.

### `No analyses are enabled`

You need at least one configured section under `plugins:`.

### `plot-all` runs but expected figures are missing

Check that the corresponding comparison JSON files already exist under
`comparison/<analysis>/result.json` and use
`polyzymd compare plot-all --list-available` to verify the enabled plot types.

### `polyzymd compare run` fails for an experimental metric

Run the prerequisite analysis first. For example, `binding_free_energy` and
`polymer_affinity` depend on cached contact-derived data, so you usually run:

```bash
polyzymd compare run contacts
polyzymd compare run binding_free_energy
```

## See Also

- [Tutorial: Analyze a Study from Finished Simulations](../tutorials/analysis_complete_workflow.md)
- [Comparison and Plotting Reference](../reference/analysis_comparison_reference.md)
- [Statistical Best Practices for Analysis](../explanation/analysis_statistics_best_practices.md)
