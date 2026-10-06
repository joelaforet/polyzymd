# Quickstart example: a protein in water

Trp-cage (PDB 1L2Y, model 1), cleaned with `polyzymd clean-pdb`, in TIP3P
water with 0.15 M NaCl. The run is a few picoseconds, so it finishes on a
laptop CPU in about two minutes.

| File | What it is |
|---|---|
| `trpcage.pdb` | The protein, chain A, with hydrogens |
| `config.yaml` | OpenMM, CPU platform |
| `config_gromacs.yaml` | The same system on GROMACS |

```bash
polyzymd validate -c config.yaml
polyzymd run -c config.yaml -r 1
polyzymd study init study --condition "Water=config.yaml" --equilibration 0ns
polyzymd analyze rg --study study
```

For a real study, lengthen `duration` (in ns) and run on a GPU. The test
`tests/test_quickstart_example.py` runs these commands for both engines.
