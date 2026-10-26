# Co-solvent example: a protein in 2 M urea

- **Shows:** Trp-cage (from `../quickstart/trpcage.pdb`) in 2 M urea with
  TIP3P water and 0.15 M NaCl, on OpenMM.
- **You learn:** how `solvent.co_solvents` adds a library co-solvent by name,
  with its shipped charges (`charges: {method: preset}`), and how
  `concentration` sets the number of molecules.
- **Time:** about 1.5 minutes on a laptop CPU, build included (6,738 atoms).
- **Look at** `build_provenance.md` in the replicate folder:
  `composition.counts` in the physics audit (75 urea, 1941 waters) and
  `composition.density` (about 1.12 g/mL: a co-solvent given as a
  concentration is added on top of the water).

| File | What it is |
|---|---|
| `config.yaml` | OpenMM, CPU platform |

The tutorial "Set up your own system" adds the same urea condition (step 5);
the configuration reference ("Co-solvents") gives every key and the formula
for the number of molecules. From the `examples/` folder, in the `build`
environment:

```bash
polyzymd project init paper --study trpcage
polyzymd study add-condition "Urea 2 M" --config cosolvent/config.yaml --study paper/trpcage
polyzymd validate -c paper/trpcage/conditions/urea_2_m/config.yaml
polyzymd run -c paper/trpcage/conditions/urea_2_m/config.yaml -r 1
# replace "analyses: {}" in paper/project.yaml with:
#   analyses:
#     rg: {}
polyzymd analyze --project paper
```

The run goes into `paper/runs/trpcage/urea_2_m/`. To compare with water, add
`../quickstart/config.yaml` to the same study first, as the `Water`
condition: the first condition of a study is its control.

To make it a real run, lengthen production `duration` (ns) and the
equilibration, and run on a GPU (`openmm.platform: CUDA`; see the how-to
"Run OpenMM on other hardware"). For another co-solvent, give a library name
(such as `dmso`) or a `smiles` with a charge method other than `preset`.
