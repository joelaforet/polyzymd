# Glossary

The docs use each of these terms with one meaning only.

## Studies and simulations

```{glossary}
project
  One paper. A project folder holds one study for each protein or system in
  the paper. Make one with `polyzymd project init`. See
  {doc}`../explanation/projects`.

study
  One protein or system, simulated under several conditions. A study folder
  holds the conditions, `study.yaml`, the stored results and the figures. Make
  one with `polyzymd study init`. See {doc}`../explanation/study_folders`.

condition
  One simulated system: one `config.yaml`. Two conditions of a study differ
  in a polymer, a co-solvent, a temperature or another setting.

replicate
  One independent simulation of a condition. The replicate number seeds the
  starting structure (the PACKMOL placement of molecules and, with polymers,
  the random draw of each chain's monomer sequence), the initial velocities
  and the thermostat noise, on both engines.

simulation folder
  The folder that `polyzymd init` makes: a template `config.yaml` and a
  `structures/` folder. It is a place to write one config. It is not a
  {term}`project`.

replicate folder
  The folder of one replicate's output, named by `output.naming_template`
  (for example `trpcage_300K_run1/`). It holds the built system and one
  `production_N/` folder for each production {term}`segment`.

segment
  One part of a production simulation. A long simulation on a cluster runs as
  a chain of SLURM jobs, one segment each. Segment N writes `production_N/`.

run
  One named entry of `analyses:` in `study.yaml`. Its stored results go to
  `results/<run>/`. Two runs can use the same analysis with different
  settings.

result
  One quantity that an analysis reports, such as `mean_rg`. An analysis can
  report several results. `polyzymd analyze --run NAME` selects the result to
  report.

config hash
  The first 16 hexadecimal characters of the SHA-256 of the config fields
  that decide what was simulated. These fields are the enzyme and its PDB
  content, the temperature and pressure, the naming template, the substrate,
  the polymers and the co-solvents. Paths, simulation phases and force
  fields are not in the hash.
```

## Force fields and system building

```{glossary}
OpenFF
  The Open Force Field toolkit. PolyzyMD uses it to read molecules and to
  assign force-field parameters. See the
  [OpenFF toolkit docs](https://docs.openforcefield.org/projects/toolkit/).

Interchange
  The OpenFF object that holds a parametrized system: topology, force-field
  parameters, positions and box. PolyzyMD exports it to OpenMM or GROMACS
  files. See the
  [Interchange docs](https://docs.openforcefield.org/projects/interchange/).

AM1-BCC
  A partial-charge method: a semi-empirical AM1 calculation with bond charge
  corrections. In PolyzyMD it needs AmberTools, which the default
  environments do not include.

NAGL
  The OpenFF graph neural network that predicts AM1-BCC partial charges. It
  needs no quantum-chemistry program. It is the default `charge_method`.
  See the [NAGL docs](https://docs.openforcefield.org/projects/nagl/).

PACKMOL
  The program that places polymer, co-solvent and solvent molecules around
  the protein. Its `tolerance` is the smallest distance, in Å, between atoms
  of two different molecules.

PME
  Particle mesh Ewald. A method that computes the long-range electrostatic
  energy in a periodic box.
```

## Statistics

```{glossary}
equilibration window
  The time at the start of each production simulation that an analysis
  discards. Set it with `--eq` or `equilibration:` in `study.yaml`.

statistical inefficiency
  The factor g by which correlation between frames reduces the information in
  a time series. Uncorrelated frames give g = 1.

n_eff
  The effective number of independent samples in a time series: the number
  of frames divided by the {term}`statistical inefficiency` g.

95 % confidence interval
  The interval that PolyzyMD reports for a condition mean. It uses Student's
  t distribution over the replicate means, so n is the number of replicates.

Benjamini-Hochberg
  A correction of p values for many tests. It controls the expected fraction
  of false discoveries among the tests that it reports as significant.
  PolyzyMD applies it within each family of tests and reports `p_adj`.
```
