# Run a simulation longer

Use this page to continue a finished OpenMM replicate for more nanoseconds.
The longer run is the same simulation: it continues from the saved state of
the last segment, with the same build, thermostat and frame spacing. Its new
frames join the old ones as one trajectory.

This works for OpenMM runs only. PolyzyMD cannot run a GROMACS replicate
longer. For GROMACS, run new replicates with the longer `duration`. You can
also extend a run yourself with `gmx convert-tpr -extend` and
`gmx mdrun -cpi`, but PolyzyMD does not record or check such a run.

An interrupted run, which stopped before the end of its production,
continues on both engines when you run the same command again. See
[Stop and continue a local run](hardware_platforms.md#stop-and-continue-a-local-run)
for OpenMM and
[Stop and continue a local run](run_gromacs.md#stop-and-continue-a-local-run)
for GROMACS.

## Change the config

The config belongs to the condition, not to one replicate. Raising its
length makes every replicate of the condition, and every replicate you run
later, run to the new length.

1. Open the condition's `config.yaml`.
2. Look at `output.naming_template`. Configs that `polyzymd study
   add-condition --new` writes leave the length out of the folder name; skip
   this step for them. A template with `{duration}`, such as the default of
   a config that sets no template,
   `{enzyme}_{substrate}_{polymer_type}_{duration}ns_{temperature}K_run{replicate}`,
   puts the production length in the replicate's folder name, and a new
   `duration` would point PolyzyMD at a new, empty folder. Write the current
   length into the template in place of `{duration}`, for example
   `{enzyme}_{substrate}_{polymer_type}_1000ns_{temperature}K_run{replicate}`.
   The replicate keeps its folder, and PolyzyMD accepts this change.
3. Raise `simulation_phases.production.duration` to the new total, for
   example from `1000.0` to `2000.0` (ns).
4. Raise `simulation_phases.production.samples` by the same factor, for
   example from `2500` to `5000`. The frame spacing is
   `duration / samples`, so it stays the same. A run that would change it is
   refused, and the message gives the `samples` value that keeps it.
5. Change nothing else.

## Continue on a workstation

Run the same command as before:

```bash
pixi run -e sim-cuda-12-6 polyzymd run -c config.yaml -r 1
```

`polyzymd run` prints `The config changes only the production length` and
runs a new segment from the end of the last one to the new total. Any other
change to the config is refused.

To preview what the run would do, add `--dry-run`. It says, per replicate,
whether it would continue from step N to M, is complete, or would be refused
and why. `polyzymd status` shows a replicate that finished its old length as
`EXTENDABLE`.

## Continue on a SLURM cluster

Submit the replicate again:

```bash
polyzymd submit -c config.yaml -r 1 --preset <preset> --pixi-env sim-cuda-12-4
```

Each job continues from the last segment and resubmits itself until the new
total is reached. `submit` skips complete replicates, and refuses one whose
next segment would change the frame spacing. `submit --dry-run` says what
each replicate's chain would do.

A replicate started on a workstation can continue on a cluster, and the other
way round. Copy its whole folder to the place the config's
`scratch_directory` names there. Each segment continues from the last
segment's portable state file (`production_N_state.xml`, or
`interrupted_state.xml`), which any OpenMM version reads. Only a segment that
left nothing but a binary checkpoint needs the same OpenMM build.

On a cluster, `submit` does not compare a started replicate with its build.
Each new segment takes its thermostat, friction and barostat from the
parameters file of the segment before it, not from the config. So a change to
another setting is not refused there: it is ignored, or it makes the chain
inconsistent. Edit only the template, `duration` and `samples`.

## Check the result

```bash
polyzymd status -c config.yaml
polyzymd analyze rg -c config.yaml
```

`status` shows the new total. The analysis reads every segment as one
trajectory, with frame times that continue from the old ones.

Stored analysis results are measured again. Each stored replicate result
records the trajectory files it read, with their sizes and hashes, and a run
that has a new segment no longer matches it. Changing the folder template
also changes the condition's config hash, which has the same effect.

To compare a longer replicate with replicates of the old length over the same
time, give `--until` (for example `--until 1000ns`), or set `until:` in the
study file.

## Limits

- Changing the temperature, the ensemble or another setting part-way is not
  supported. To study such a change, run it as a new condition.
- Starting a new condition from the final state of a finished replicate is
  not supported yet. Build and equilibrate the new condition from its own
  config instead (`polyzymd study add-condition NAME --from <condition>`).
  To start from that state anyway, run your own protocol outside PolyzyMD
  on the finished replicate's `production_N_state.xml`, and analyse the
  result as {doc}`own_protocol` describes.
