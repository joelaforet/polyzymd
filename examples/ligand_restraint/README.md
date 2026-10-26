# Ligand and restraint example: T4 lysozyme L99A with benzene

- **Shows:** T4 lysozyme L99A (PDB 181L) with benzene at its crystal pose in
  the buried cavity, in TIP3P water with 0.15 M NaCl, on OpenMM. A
  flat-bottom distance restraint keeps the benzene in the cavity, and a first
  equilibration stage holds the protein and benzene heavy atoms in place.
- **You learn:** how a `small_molecule` solute reads an SDF, how a restraint
  names its atoms by `component` and selection (the benzene carbon is `C1x`,
  not `C1`), and how a PDB entry is prepared for OpenFF.
- **Time:** about 3.5 minutes on a laptop CPU, build included (33,864
  atoms). With a laptop GPU, about 1.5 minutes: 1 minute to build, 30 s to
  run (see the end of this file).
- **Look at** `build_provenance.md` in the replicate folder: the `benzene`
  row of Components (NAGL charges), and `restraint.atom_identity` and
  `restraint.parameters` in the physics audit, which name the two restrained
  atoms (`ALA99:CB`, `BNZ1:C1x`). The run log prints
  `Applied flat_bottom restraint 'benzene_in_cavity'`.

| File | What it is |
|---|---|
| `config.yaml` | OpenMM, CPU platform, with the restraint |
| `structures/t4l_clean.pdb` | The protein, prepared from 181L (2603 atoms, net charge +8) |
| `structures/benzene.sdf` | The benzene at its crystal coordinates, with hydrogens (12 atoms) |
| `add_terminal_oxt.py`, `ligand_to_sdf.py` | The two scripts that made them |

The tutorial "Prepare a PDB for OpenFF and PolyzyMD" (section "A chain that
ends early, and a crystal ligand") explains the preparation, and the how-to
"Add distance restraints" the restraint. From the `examples/` folder, in the
`build` environment:

```bash
polyzymd project init paper --study t4l
polyzymd study add-condition Benzene --config ligand_restraint/config.yaml --study paper/t4l
polyzymd validate -c paper/t4l/conditions/benzene/config.yaml
polyzymd run -c paper/t4l/conditions/benzene/config.yaml -r 1
# replace "analyses: {}" in paper/project.yaml with:
#   analyses:
#     rmsd: {}
polyzymd analyze --project paper
```

The run goes into `paper/runs/t4l/benzene/`.

## How the input files were made

From RCSB entry 181L (downloaded on 2026-10-25, SHA-256
`77018feaaa65bb22dea47c784e8c059b0ccc09cd6dc7442b79cce83f3170985f`), with
PolyzyMD 1.3.0, OpenMM 8.4, PDBFixer and RDKit of the `build` environment,
following the tutorial exactly:

```bash
curl -L https://files.rcsb.org/download/181L.pdb -o structures/181L.pdb
grep -v "^HETATM" structures/181L.pdb > structures/181L_protein.pdb
python add_terminal_oxt.py structures/181L_protein.pdb structures/181L_oxt.pdb
polyzymd clean-pdb -i structures/181L_oxt.pdb -o structures/t4l_clean.pdb --ph 7.0
python ligand_to_sdf.py structures/181L.pdb BNZ "c1ccccc1" structures/benzene.sdf
```

`grep` removes the waters, ions, additives and the benzene. 181L has no
coordinates for residues 163 and 164, so `add_terminal_oxt.py` adds the
terminal oxygen of Lys 162 (and nothing else). `clean-pdb` adds hydrogens at
pH 7. `ligand_to_sdf.py` takes the bond orders from the SMILES and keeps the
crystal coordinates of benzene `BNZ A 400`. OpenFF reads the result as one
molecule of 2603 atoms with net charge +8.

The restraint atoms: Ala 99 CB lines the cavity the L99A mutation makes and
is 3.8 Å from the first benzene carbon in the crystal (C1, `C1x` in the built
system). The flat bottom is 6.5 Å, so the restraint acts only if the benzene
starts to leave the cavity.

To run on a local NVIDIA GPU, set `platform: CUDA` (and `precision: mixed`)
under `openmm:` in the condition's `config.yaml`, build in the `build`
environment and run in a `sim-cuda-*` environment that your driver supports
(see the how-to "Run OpenMM on other hardware"):

```bash
pixi run -e build polyzymd build -c paper/t4l/conditions/benzene/config.yaml -r 1
pixi run -e sim-cuda-12-6 polyzymd run -c paper/t4l/conditions/benzene/config.yaml -r 1
```

To make it a real run, lengthen production `duration` (ns) and the
equilibration, and run three or more replicates on a GPU. Distance restraints
need `engine: openmm`.
