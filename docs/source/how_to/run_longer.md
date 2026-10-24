# Run a simulation longer

Use this page to continue a finished OpenMM replicate for more nanoseconds.
The longer run is the same simulation: it continues from the saved state of
the last segment, with the same build, thermostat and frame spacing. Its new
frames join the old ones as one trajectory.

## Raise the production length

1. Open the replicate's `config.yaml`.
2. Raise `simulation_phases.production.duration` to the new total, for
   example from `1000.0` to `2000.0` (ns).
3. Raise `simulation_phases.production.samples` by the same factor, for
   example from `2500` to `5000`. The frame spacing is
   `duration / samples`, so it stays the same. A run that would change it is
   refused, and the message gives the `samples` value that keeps it.
4. Change nothing else. A started production continues only with the config
   it was built with, apart from these two settings.

## Continue on a workstation

Run the same command as before:

```bash
pixi run -e sim-cuda-12-6 polyzymd run -c config.yaml -r 1
```

`polyzymd run` prints `The config changes only the production length` and
runs a new segment from the end of the last one to the new total.

## Continue on a SLURM cluster

Submit the replicate again:

```bash
polyzymd submit -c config.yaml -r 1 --preset <preset>
```

Each job runs `polyzymd run-segment`, which continues from the last segment.
`polyzymd check-progress` reads the new total from the config, so the chain
resubmits itself until it reaches it.

On a cluster, `submit` does not compare a started replicate with its build.
Each new segment takes its thermostat, friction and barostat from the
parameters file of the segment before it, not from the config. So a change to
another setting is not refused there: it is ignored, or it makes the chain
inconsistent. Edit only `duration` and `samples`.

## Check the result

```bash
polyzymd check-progress -c config.yaml -r 1
polyzymd analyze rg -c config.yaml
```

`check-progress` reports the new total and one more segment. The analysis
reads every segment as one trajectory, with frame times that continue from
the old ones.

## Limits

- Only OpenMM runs can be continued this way.
- Changing the temperature, the ensemble or another setting part-way is not
  supported. To study such a change, run it as a new condition.
