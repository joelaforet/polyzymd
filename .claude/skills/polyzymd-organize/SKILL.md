---
name: polyzymd-organize
description: Organize PolyzyMD simulations into a project, studies and conditions, and keep every file where PolyzyMD expects it. Use when a user starts a new set of simulations, adds a condition (a new polymer, temperature, solvent or ligand), wants to compare runs that already exist, or asks where configs, inputs, runs, logs or results belong. Use the polyzymd project and study commands; never make condition folders by hand.
---

# polyzymd-organize: one project, studies, conditions

PolyzyMD keeps a paper's simulations in one **project** folder. A project holds
**studies**; a study holds **conditions**; a condition is one `config.yaml`
with its input files, run as several replicates. Build and run a condition
with the `polyzymd-simulate` skill; analyse a study with `polyzymd-analyze`.

## 1. Where every file lives

```
my_paper/                         project (git repository)
├── project.yaml                  its studies, shared analyses, publishing metadata
├── analyses/  stats/  figures/   project-wide analysis code, statistics, figures
├── <study>/                      one study
│   ├── study.yaml                conditions (control first), structures, regions, window
│   ├── data.local.yaml           where this machine keeps each condition's runs (git-ignored)
│   ├── structures/               reference structures the analyses name
│   ├── conditions/<condition>/
│   │   ├── config.yaml           the condition; comments say why it differs
│   │   └── structures/           every input file the config names, copied in
│   └── results/                  stored analysis results
└── runs/<study>/<condition>/     default runs, job scripts and logs (git-ignored)
```

- Configs, input files, `study.yaml`, `project.yaml`, analysis code and stored
  results are tracked by git. Trajectories are not.
- A config's runs go to `output.scratch_directory` when it is set, otherwise
  to `runs/<study>/<condition>/`. Job scripts and SLURM logs go to
  `output.projects_directory`.
- Never keep trajectories in the project, and never commit `runs/`.

## 2. Decide the studies first

A study is one analysis frame: one residue numbering, one set of reference
structures and regions, one equilibration window and **one control** (its
first condition). Every analysis compares each condition with that control.
Choose the control from the user's question, then:

- A different protein: a new study.
- The question is "what does X do?" (polymer, co-solvent, ligand) at one
  temperature: one study per temperature, with that temperature's no-X run as
  the control. Example: `lipa363` and `lipa500`, each with its own No Polymer.
- The question is "what does temperature do?": one study whose conditions are
  the temperatures, with the reference temperature as the control.
- The same runs answer two questions: add the condition to both studies with
  `--config`; the copy reuses the runs (`docs/source/how_to/reuse_condition.md`).

Ask the user when the control is not obvious. Name studies and conditions as
the user would say them (`"SBMA-EGMA 100:0"`); the folder name is derived.

## 3. Start a project or add to one

```bash
polyzymd project init my_paper --study lipa363 --study lipa500   # skeleton, git repository
cd my_paper
polyzymd study add-condition "No Polymer" --new --study lipa500            # template to fill in
polyzymd study add-condition "No Polymer" --config ../old/config.yaml --study lipa500  # copy a config
polyzymd study add-condition "SBMA-EGMA 100:0" --from "No Polymer" --study lipa500    # copy, then edit
polyzymd project add-study calb343 --project .
```

`--config` and `--from` copy every input file the config names into the
condition's `structures/` and rewrite the paths. Then, after each change:

1. Read the copied `config.yaml`. Check `output.scratch_directory` (step 4)
   and that the comments explaining the condition are still there.
2. Record what varies in `study.yaml`:
   `"SBMA-EGMA 100:0": {config: conditions/sbma-egma_100_0/config.yaml, factors: {sbma_fraction: 1.0}}`.
3. `polyzymd validate -c <study>/conditions/<condition>/config.yaml`
4. `polyzymd study check <study>`: control, replicates found, uncommitted inputs.
5. Commit: `git add -A <study> && git commit -m "<study>: add <condition>"`.

## 4. Runs on a cluster: set scratch_directory

Trajectories are large. On a cluster, every condition's config must set
`output.scratch_directory` to scratch storage before you build or submit.

- `add-condition --config` keeps `scratch_directory` only when the copied
  config already has runs there (it says "this config reuses them"). For a
  config with no runs yet, the copy's `scratch_directory` becomes null and the
  runs would go into `runs/` inside the project. Set it back, then validate.
- `--new` and `--from` also leave it null: set it.
- After moving or downloading runs, `polyzymd study locate DIR` records their
  folder in `data.local.yaml`; do not edit configs to point at them.

## 5. Bring in runs that already exist

Make a new project next to the old folders; do not move or rename anything a
running job uses.

1. `polyzymd project init` with one study per old study or enzyme.
2. `polyzymd study add-condition LABEL --config OLD/config.yaml --study S` for
   each condition, control first. Configs with runs keep their
   `scratch_directory`, so `study check` finds the replicates.
3. Running SLURM chains keep reading their old config. Do not edit the old
   config or the copy while a chain runs; `study check` warns when the two
   differ. Submit new replicates from the project's copy.
4. Follow `docs/source/how_to/move_studies_into_project.md` for structures,
   regions and stored results.

## 6. Monitor

```bash
polyzymd status --format agent --all my_paper --unfinished   # every condition of the project
polyzymd study check my_paper/lipa500
polyzymd project check my_paper
```

## 7. Rules

- Make projects, studies and conditions only with `polyzymd project` and
  `polyzymd study`; never copy condition folders by hand.
- One condition, one config. Write why it differs from the control in the
  config's top comment.
- Never edit the config of a condition whose chains are running, and never
  change a condition's physics after its first replicate; add a new condition.
- Keep inputs inside the condition (`structures/`); never point a config at a
  file outside the project.
- When the commands get in the way (a path rewritten, a warning that does not
  match what happened, a step the CLI cannot do), stop and tell the user what
  happened and what you did instead.
