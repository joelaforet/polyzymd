# Run a simulation longer

Use this page to continue a finished OpenMM replicate for more nanoseconds.
The longer run is the same simulation: it continues from the saved state of
the last segment, with the same build, thermostat and frame spacing. Its new
frames join the old ones as one trajectory.

This works for OpenMM runs only. PolyzyMD cannot run a GROMACS replicate
longer. For GROMACS, run new replicates with the longer `duration`. You can
also extend a run yourself with `gmx convert-tpr -extend` and
`gmx mdrun -cpi`, but PolyzyMD does not record or check such a run.

## Change the config

1. Open the replicate's `config.yaml`.
2. Look at `output.naming_template`. The default,
   `{enzyme}_{substrate}_{polymer_type}_{duration}ns_{temperature}K_run{replicate}`,
   puts the production length in the replicate's folder name. A new
   `duration` would then point PolyzyMD at a new, empty folder. Write the
   current length into the template in place of `{duration}`, for example
   `{enzyme}_{substrate}_{polymer_type}_1000ns_{temperature}K_run{replicate}`.
   The replicate keeps its folder, and PolyzyMD accepts this change. Skip
   this step if the template has no `{duration}`.
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

## Continue on a SLURM cluster

Submit the replicate again:

```bash
polyzymd submit -c config.yaml -r 1 --preset <preset>
```

Each job continues from the last segment and resubmits itself until the new
total is reached.

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
