# Set Up Equilibration Stages

Use this guide to configure PolyzyMD's required staged equilibration workflow
with heating, restrained relaxation, or a final unrestrained equilibration
step.

## When to use staged equilibration

Multi-stage equilibration is useful when:

- the initial packed system has bad contacts
- you want gradual heating instead of starting directly at the target temperature
- polymers or ligands should relax while the protein stays restrained
- you want a final free equilibration stage before production

## Basic rule

PolyzyMD now requires `simulation_phases.equilibration_stages`.

Even if your protocol is minimal, represent it as one or more named stages.

## Simple three-stage example

```yaml
simulation_phases:
  equilibration_stages:
    - name: "heating"
      samples: 50
      ensemble: "NVT"
      temperature_start: 60.0
      temperature_end: 293.0
      temperature_increment: 1.0
      temperature_interval_steps: 600  # duration is derived
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0
        - group: "polymer_heavy"
          force_constant: 4184.0

    - name: "polymer_relaxation"
      duration: 1.0
      samples: 100
      ensemble: "NVT"
      temperature: 293.0
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0

    - name: "free_equilibration"
      duration: 0.5
      samples: 100
      ensemble: "NPT"
      temperature: 293.0

  production:
    ensemble: "NPT"
    duration: 100.0
    samples: 2500
    time_step: 2.0
    thermostat: "LangevinMiddle"
    thermostat_timescale: 1.0
    barostat: "MC"
    barostat_frequency: 25
```

## How to choose stage settings

### Heating stage

Use a heating stage when you want a gentler start for a crowded or fragile
system. Set `temperature_start`, `temperature_end`, and
`temperature_increment` in K, plus `temperature_interval_steps` in MD steps. Do
not set `duration`; PolyzyMD calculates the number of updates needed to reach
the endpoint and derives the duration. If the final update would overshoot,
PolyzyMD shortens only that final temperature increase. Validation and
simulation logs report the derived duration before the stage runs.

The ramp ends when the target reaches `temperature_end`; it does not include a
hold at that temperature. Add a following constant-temperature stage when the
system should equilibrate at the endpoint. OpenMM updates the thermostat target
at each requested step interval, while GROMACS interpolates between the same
temperature update points over the same derived duration.

### Restrained relaxation stage

Use `position_restraints` when one part of the system should stay close to its
reference coordinates while another part relaxes.

Available groups include:

- `protein_heavy`
- `protein_backbone`
- `protein_calpha`
- `ligand_heavy`
- `polymer_heavy`
- `solvent`, `water_only`, `ions_only`, `cosolvents_only`

### Free equilibration stage

Omit `position_restraints` when you want the stage to run fully unrestrained.

## Common recipes

### Enzyme-only system

```yaml
simulation_phases:
  equilibration_stages:
    - name: "heating"
      ensemble: "NVT"
      temperature_start: 60.0
      temperature_end: 300.0
      temperature_increment: 1.0
      temperature_interval_steps: 600
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0

    - name: "free_equilibration"
      duration: 1.0
      ensemble: "NPT"
      temperature: 300.0
```

### Gradual restraint release

```yaml
simulation_phases:
  equilibration_stages:
    - name: "heating"
      ensemble: "NVT"
      temperature_start: 60.0
      temperature_end: 300.0
      temperature_increment: 1.0
      temperature_interval_steps: 600
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0

    - name: "relax_1"
      duration: 0.5
      ensemble: "NVT"
      temperature: 300.0
      position_restraints:
        - group: "protein_heavy"
          force_constant: 2092.0

    - name: "relax_2"
      duration: 0.5
      ensemble: "NVT"
      temperature: 300.0
      position_restraints:
        - group: "protein_backbone"
          force_constant: 1046.0

    - name: "free_equilibration"
      duration: 0.5
      ensemble: "NPT"
      temperature: 300.0
```

## Validate the setup

Run:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd build -c config.yaml --dry-run
```

## Outputs to expect

Each stage writes its own trajectory, state-data CSV, topology PDB, and
checkpoint files. The stage name becomes part of the output filename, which
makes it easier to inspect temperature ramps or restraint-release steps later.

## Troubleshooting

### simulation blows up during heating

Try a smaller `time_step`, a lower starting temperature, or stronger initial
position restraints.

### temperature overshoots the target

Use a smaller `temperature_increment` or a larger `temperature_interval_steps`.

### protein moves too much while polymers relax

Keep `protein_heavy` restraints for longer or reduce them more gradually.

## Related pages

- first-run tutorial: {doc}`../get_started/quickstart`
- distance restraints: {doc}`restraints`
- config schema: {doc}`../reference/configuration`

<!-- IMAGE OPPORTUNITY: Add a three-panel schematic for `heating -> restrained
relaxation -> free equilibration`, with example restraint groups called out in
each panel. -->
