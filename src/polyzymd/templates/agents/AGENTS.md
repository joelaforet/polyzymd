# Working in this PolyzyMD project

This folder holds molecular dynamics simulations organized by PolyzyMD. Use the
`polyzymd` commands below for every step; do not write your own job scripts,
trajectory loaders, or loops over folders.

## What the folders mean

A project is one paper: `project.yaml` lists its studies. A study is a set of
conditions compared with each other (one protein, one residue numbering). A
condition is one simulation setup, a `config.yaml`. A replicate is one
independent run of a condition. A lone study has the same layout without
`project.yaml`.

| Path | What it holds |
|------|---------------|
| `project.yaml` | The studies, the analyses every study runs, publishing metadata |
| `analyses/` | Measurement functions, named as `file.py:function` |
| `stats/`, `figures/` | Statistics scripts and figure scripts that read stored results |
| `<study>/study.yaml` | Conditions, structures, regions, equilibration window, analyses |
| `<study>/conditions/<condition>/config.yaml` | One condition's simulation setup |
| `<study>/results/<run>/` | Reports and figures written by `polyzymd analyze` |
| `<study>/sandbox/` | One-off checks and exploratory scripts |
| `runs/<study>/<condition>/` | Trajectories of each replicate; never committed |

## The commands

```bash
polyzymd project check .                          # every study: conditions, runs, problems
polyzymd study check <study>                      # one study
polyzymd study add-condition "LABEL" --from OTHER --study <study>   # or --config FILE, or --new
polyzymd submit -c <study>/conditions/<condition>/config.yaml -r 1-5 --preset <cluster>
polyzymd build -c <study>/conditions/<condition>/config.yaml         # build without running
polyzymd status --all . --format agent            # progress of every replicate
polyzymd analyze --list                           # the analyses PolyzyMD ships
polyzymd analyze --study <study> --format agent   # every analysis study.yaml lists
polyzymd analyze --project . --format agent       # every study of the project
polyzymd analyze <run> --study <study>            # one run that study.yaml lists
polyzymd project freeze .                         # freeze for publication (study freeze for one)
polyzymd agent sync                               # refresh these agent files
polyzymd project list --format agent              # your other projects and where they are
```

## Rules

- Run `polyzymd analyze --list` before writing any MDAnalysis code. Most
  questions are answered by a shipped analysis.
- To measure something new on every replicate, write a function of an
  MDAnalysis `Universe` in `analyses/` and list it in `study.yaml` (or
  `project.yaml`) under `analyses:` with `function: analyses/file.py:name`,
  then run `polyzymd analyze --study <study>`. In Python, use
  `Study.timeseries` instead.
- Never loop over replicates or study folders by hand, and never build file
  paths to trajectories yourself. PolyzyMD finds every replicate, applies the
  equilibration window, and records what it read.
- Put one-off scripts only in `<study>/sandbox/`, never in the project root
  or the current folder. A script that becomes a real analysis moves to
  `analyses/`.
- Never choose a force field, charge model, water model, or temperature for
  the user. Ask, and copy the values of an existing condition when told to.
- Commit inputs (configs, structures, `study.yaml`, `project.yaml`, analysis
  functions) after each decision, one commit per decision, with a message
  that says why.
- PolyzyMD never uploads anything. Uploading to Zenodo or another archive is
  the author's step, done by hand after `freeze`.
- If this project lives on a network mount (sshfs), read and edit single
  files through the mount, but run git, `polyzymd`, and anything that walks
  the folder tree on the machine that hosts the files. Never run
  `git status`, a recursive `grep`, or a glob over the mount.
- Commands print warnings and problems; act on them before going on. The
  full log of each command is in `logs/`.
