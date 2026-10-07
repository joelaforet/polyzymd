# Quickstart example: a protein in water

Trp-cage (PDB 1L2Y, model 1), cleaned with `polyzymd clean-pdb`, in TIP3P
water with 0.15 M NaCl. The run is a few picoseconds, so it finishes on a
laptop CPU in about two minutes.

| File | What it is |
|---|---|
| `trpcage.pdb` | The protein, chain A, with hydrogens |
| `config.yaml` | OpenMM, CPU platform |
| `config_gromacs.yaml` | The same system on GROMACS |

The quickstart tutorial adds this config to a new project as a condition,
then runs and analyzes it. From the folder that holds this one:

```bash
polyzymd project init paper --study trpcage
polyzymd study add-condition Water --config quickstart/config.yaml --study paper/trpcage
polyzymd validate -c paper/trpcage/conditions/water/config.yaml
polyzymd run -c paper/trpcage/conditions/water/config.yaml -r 1
# list rg under analyses: in paper/project.yaml, then:
polyzymd analyze --project paper
```

The run goes into `paper/runs/trpcage/water/`, which git ignores.

For a real study, lengthen `duration` (in ns) and run on a GPU. The test
`tests/test_quickstart_example.py` runs these commands for both engines.
