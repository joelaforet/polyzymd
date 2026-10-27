# Polymer example: a protein with polymer chains built from SMILES

- **Shows:** Trp-cage (from `../quickstart/trpcage.pdb`) with two sulfobetaine
  methacrylate (SBMA) 5-mers, built from star-sited SMILES by mBuild and
  charged with NAGL, in TIP3P water with 0.15 M NaCl. One config per engine.
- **You learn:** how a `type: polymer` entry with `generation_mode: dynamic`
  turns a monomer SMILES into chains, and how the same system runs on OpenMM
  and on GROMACS.
- **Time:** on a laptop CPU, build included: OpenMM about 1.5 minutes,
  GROMACS (2 threads) under 1 minute. The system has 11,244 atoms; making
  the chain takes a few seconds of the build.
- **Look at** `build_provenance.md` in the replicate folder: the `polymer`
  row of Components (NAGL model), `composition.polymer_sequences` and
  `geometry.polymer_overlap` in the physics audit, and the mBuild version and
  commit under "Inputs and seeds".

| File | What it is |
|---|---|
| `config.yaml` | OpenMM, CPU platform |
| `config_gromacs.yaml` | The same system on GROMACS (needs `gmx` on `PATH`) |

The chemistry and settings are those of the tutorial "Compare a polymer
condition with water"; the how-to "Generate polymers from SMILES" explains
every key. From the `examples/` folder, in the `build` environment:

```bash
polyzymd project init paper --study trpcage
polyzymd study add-condition SBMA --config polymer/config.yaml --study paper/trpcage
polyzymd validate -c paper/trpcage/conditions/sbma/config.yaml
polyzymd run -c paper/trpcage/conditions/sbma/config.yaml -r 1
# replace "analyses: {}" in paper/project.yaml with:
#   analyses:
#     rg: {}
#     contacts: {}
polyzymd analyze --project paper
```

For GROMACS, give `--config polymer/config_gromacs.yaml` instead. The GROMACS
run notes that barostat `MC` runs as c-rescale: GROMACS has no Monte Carlo
barostat (see the GROMACS and OpenMM reference page). `add-condition` copies `trpcage.pdb` into the condition, and the
build writes the chain into the condition's `.polymer_cache/`. The run goes
into `paper/runs/trpcage/sbma/`.

To make it a real run, lengthen production `duration` (ns) and the
equilibration, remove `packing.padding: 0.5` (the default places the chains
2.0 nm from the protein), raise `count`, and run on a GPU
(`openmm.platform: CUDA`; see the how-to "Run OpenMM on other hardware").
