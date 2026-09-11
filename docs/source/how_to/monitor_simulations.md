# Monitor a Simulation Campaign

This guide shows how to answer three questions about a set of running
PolyzyMD chains with one command: which replicates have finished, which have
died and will not restart on their own, and how long the rest will take. It
also shows how to hand that command to an AI coding agent as a reusable skill.

## Prerequisites

- Simulations submitted with `polyzymd submit` (see {doc}`hpc_slurm`), so each
  replicate has a `progress.json` in scratch and SLURM logs in
  `<projects_dir>/slurm_logs/`.
- A shell on the cluster login node with `squeue` available. On clusters with
  several schedulers, load the one your jobs run under first (on CU Boulder
  Blanca this is `ml slurm/blanca`).
- A PolyzyMD pixi environment that has the CLI, for example `build` or
  `analysis`.

## Check one system

```bash
pixi run -e build polyzymd status --format agent -c config.yaml --preset <preset>
```

```
# polyzymd status  2026-09-11 16:28 UTC  1 system(s)  5 replicate(s): 3 running, 1 queued, 1 dead

## CALB_ResorufinButyrate_none_1000ns_343K  (config.yaml)
run1   361.6/1000ns   36%  RUNNING      job 28248421 R 5:27 bgpu-shirts3  176ns/d  eta 3.6d
run3   364.0/1000ns   36%  RUNNING      job 28248422 R 1:38 bgpu-biokem1  185ns/d  eta 3.4d
run4   354.2/1000ns   35%  QUEUED       job 28248423 PD ((Priority))  eta ?
run5   708.2/1000ns   71%  RUNNING      job 28248388 R 39:13 bgpu-shirts2  324ns/d  eta 22h
run2   112.3/1000ns   11%  DEAD         no job  last: FATAL: CUDA routing failed after 3 retries [CALB_..._run2.28228845.out]

# dead chains — resume from checkpoint with:
polyzymd submit -c config.yaml -r 2 --preset blanca-shirts
```

Read the fourth column first. `RUNNING` and `QUEUED` mean a SLURM job with
this replicate's name exists, so the chain is alive even if `progress.json`
says `interrupted` (chains roll over at the wall-time limit, and a job that
started minutes ago is normal). `DEAD` means work remains and nothing is
queued or running; the `last:` field is the line from the newest SLURM log
that explains why, and the footer gives the command that resumes it from its
checkpoint. The full verdict table is in {ref}`cli-status`.

## Check a whole campaign

Point `--all` at the directory tree that holds your config files. It finds
every `config.yaml` up to three levels deep, makes a single `squeue` call for
all of them, and prints one block per system:

```bash
pixi run -e build polyzymd status --format agent \
    --all /projects/$USER/sims/CALB \
    --all /projects/$USER/sims/RML \
    --preset blanca-shirts
```

Config files that no longer match the current schema print one
`ERROR failed to load config` line and are otherwise skipped, so old
directories in the tree do not stop the report.

The header line is the campaign summary. For "how much longer", read the
`eta` column as a range across replicates rather than a single number:
preemptable GPU nodes differ several-fold in speed, and a replicate that just
restarted has only its previous segment to measure from.

## Act on dead chains

Copy the `polyzymd submit` lines from the footer and run them from the
directory that holds the config, with the correct scheduler module loaded.
`submit` calls `sbatch` directly, so if the wrong scheduler is active the
submission fails with `invalid partition specified`. Resubmission resumes from
the last checkpoint; it never rebuilds the system.

Some `last:` lines call for something other than a plain resubmit:

| `last:` line | Meaning | Action |
|---|---|---|
| `FATAL: CUDA routing failed after 3 retries` | The chain landed on nodes whose driver is too old for the pinned CUDA environment. | Resubmit. Presets exclude the known nodes; the terminal message names any others to add with `--exclude`. |
| `CONCURRENT: Another job is already running this replicate` | A duplicate chain exited to protect the running one. | Check whether the other chain is still alive before resubmitting. |
| `Segment N failed: Particle coordinate is NaN` | The physics blew up. | Do not resubmit blindly; inspect the system. |
| `Validation error: ... polymer atom(s) lie within ... of the solute` | The build refused the packed coordinates. | Rebuild; see {doc}`broken_molecules_debugging`. |

## Script it

`--format json` emits the same data for a cron job, a dashboard, or a test:

```bash
pixi run -e build polyzymd status --format json --all /projects/$USER/sims \
    | jq -r '.systems[].replicates[] | select(.verdict=="dead") | "\(.directory) \(.last_error)"'
```

## Give the command to an AI agent as a skill

Coding agents such as Claude Code, Codex, or Cursor can answer "check on my
simulations" for you, but left alone they reconstruct the answer from
`squeue`, `sacct`, and log greps, which is slow, expensive, and easy to get
wrong. A short skill file tells the agent what the question means and which
single command answers it.

A skill is a Markdown file the agent loads when a request matches its
description. For Claude Code, create
`~/.claude/skills/polyzymd-status/SKILL.md`; other agents use an equivalent
location (`AGENTS.md`, `.cursor/rules/`, or a project instructions file).
The content is the same. Adapt the paths and preset to your cluster:

```markdown
---
name: polyzymd-status
description: Answer "check the status of the running simulations", "did anything die",
  or "how much longer" for PolyzyMD chains on the cluster with one command,
  polyzymd status --format agent. Do not start from squeue, sacct, or grep.
---

When the user asks about simulation progress, completion, deaths, restarts,
or ETA, run this first and answer from its output:

    ml slurm/blanca
    cd /projects/$USER/sims
    pixi run -e build polyzymd status --format agent --preset blanca-shirts \
        --all CALB --all RML

Read the verdict column. COMPLETED and RUNNING/QUEUED need no action.
DEAD means no job is driving the replicate; report the `last:` line
verbatim and, unless it is a NaN or a validation error, run the
`polyzymd submit` command from the footer (with the scheduler module
loaded and from the config's directory). NOT_STARTED usually means the
build failed; report it and stop.

Report: one summary line (counts and ETA range), a list of DEAD and
NOT_STARTED replicates with their reason, and what you did about them.
Do not grep slurm_logs or run sacct before running the status command.
```

Three things make a skill like this work well:

1. **Name the user's phrase**, not the command. The agent matches on the
   description, so it should contain the words the user actually says.
2. **Pin the environment.** The most common failure is running `polyzymd`
   from the wrong Python or with the wrong scheduler loaded. Put the exact
   `pixi run -e ...` and module-load lines in the skill.
3. **Say what not to do.** Agents default to exploration. An explicit
   "do not start from squeue/grep" is what saves the tokens.

Keep the skill next to the cluster-specific facts it depends on (preset
name, config root, scheduler module) and update it when those change.
