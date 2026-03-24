# Documentation Image Plan

This file tracks where images should be added after the text overhaul is stable.

## Home and Landing Pages

- `docs/source/index.rst`
  - Placement: below `Common Workflows`
  - Suggested image: one clean workflow schematic showing `Build -> Submit -> Analyze -> Compare -> Plot`
  - Purpose: help new users orient quickly without reading dense prose

- `docs/source/get_started/index.md`
  - Placement: after `Start Here`
  - Suggested image: annotated terminal screenshot of `pixi install -e build`, `pixi shell -e build`, and `polyzymd info`
  - Purpose: reassure users that the install path is short and concrete

- `docs/source/tutorials/index.md`
  - Placement: after `Featured Tutorial`
  - Suggested image: analysis pipeline figure from per-condition analysis to comparison figures
  - Purpose: make the tutorial section feel task-oriented and visual

- `docs/source/how_to/index.md`
  - Placement: between `Analysis Tasks` and `Recipe Collections`
  - Suggested image: project tree for a comparison study, with `comparison.yaml`, `results/`, and `figures/` annotated
  - Purpose: show users what artifacts to expect from the workflow

- `docs/source/reference/index.md`
  - Placement: after `API Reference`
  - Suggested image: compact module map for `config`, `builders`, `simulation`, `workflow`, `analysis`, and `compare`
  - Purpose: support scanning and improve LLM/user orientation in the reference section

- `docs/source/explanation/index.md`
  - Placement: after `Experimental Analysis Concepts`
  - Suggested image: two-panel figure showing `stable` vs `experimental` analysis families plus a chain convention diagram (`A/B/C/D+`)
  - Purpose: reinforce conceptual framing before deeper reading

- `docs/source/contributor_guide/index.md`
  - Placement: after `Extension Workflows`
  - Suggested image: high-level package architecture diagram
  - Purpose: help contributors see extension points at a glance

## Batch 2 Pages

- `docs/source/tutorials/installation.md`
  - Placement: after `Install on a local machine`
  - Suggested image: screenshot of a successful `polyzymd --help` / `polyzymd info` run inside a pixi shell
  - Purpose: show the success state clearly

- `docs/source/contributor_guide/setup.md`
  - Placement: after `What this environment is for`
  - Suggested image: a small terminal transcript showing `pytest`, `ruff`, and `make -C docs clean html`
  - Purpose: make contributor expectations obvious

- `docs/source/tutorials/analysis_complete_workflow.md`
  - Placement: after `The study we will analyze`
  - Suggested image: campaign directory tree with three conditions and a comparison workspace
  - Purpose: anchor the tutorial in a concrete example layout
  - Placement: after `Step 4: run the cross-condition comparison`
  - Suggested image: one example output figure from `polyzymd compare plot-all`
  - Purpose: give the tutorial a visual payoff before the final smoke test

- `docs/source/tutorials/analysis_compare_conditions.md`
  - Placement: after `Step 2: define a minimal comparison.yaml`
  - Suggested image: annotated `comparison.yaml` snippet or a compact infographic showing `conditions + analysis_settings + comparison_settings + plot_settings`
  - Purpose: help users parse the config without wading through raw YAML

- `docs/source/reference/analysis_comparison_reference.md`
  - Placement: near the project layout section
  - Suggested image: small reference diagram showing where result JSON files and figures are written
  - Purpose: make output locations memorable

## Batch 3 Pages

- `docs/source/tutorials/quickstart.md`
  - Placement: after `Step 1: Create a new project scaffold`
  - Suggested image: screenshot of the generated project tree from `polyzymd init`
  - Purpose: make the success state of project creation obvious for new users

- `docs/source/tutorials/hpc_slurm.md`
  - Placement: after `What the generated scripts do`
  - Suggested image: simple lifecycle figure showing `submit -> run-segment -> check-progress -> resubmit`
  - Purpose: explain the self-resubmitting job model without forcing users to read a long script first
  - Placement: near `Step 3: submit one small test job`
  - Suggested image: annotated screenshot of the top section of a generated SLURM script
  - Purpose: help users connect CLI flags to emitted `SBATCH` directives

- `docs/source/tutorials/restraints.md`
  - Placement: after `Step 4: find the right atom indices`
  - Suggested image: PyMOL screenshot of `solvated_system.pdb` with one protein atom and one ligand atom labeled
  - Purpose: show how visual atom picking maps onto YAML `selection` fields

- `docs/source/tutorials/equilibration.md`
  - Placement: after `Simple three-stage example`
  - Suggested image: three-panel schematic for `heating -> restrained relaxation -> free equilibration`
  - Purpose: help users see the role of each stage at a glance

- `docs/source/tutorials/contributing.md`
  - Placement: after `Recommended verification commands`
  - Suggested image: contributor workflow figure showing `branch -> edit -> test -> docs build -> PR`
  - Purpose: make the standard contribution loop easy to remember
  - Placement: inside `A concrete example: add a new co-solvent`
  - Suggested image: compact diagram showing `cosolvent_library.py -> generator -> solvent SDF -> dry-run build`
  - Purpose: make the example contribution path feel concrete and reproducible

- `docs/source/tutorials/architecture.md`
  - Placement: after `How data moves through the system`
  - Suggested image: left-to-right architecture diagram showing `config -> builders -> simulation/workflow -> analysis -> compare -> plots`
  - Purpose: help both contributors and advanced users orient quickly in the codebase

## Later Candidate Pages

- `docs/source/tutorials/hpc_slurm.md`
  - Suggested image: SLURM submission lifecycle or job-script anatomy

- `docs/source/tutorials/analysis_binding_preference.md`
  - Suggested image: enrichment heatmap example with stable/experimental badge

- `docs/source/tutorials/analysis_exposure_dynamics.md`
  - Suggested image: event schematic for buried -> exposed -> polymer-contact -> reburied

- `docs/source/tutorials/residue_assignment.md`
  - Suggested image: chain identity diagram for protein, substrate, polymer, and solvent
