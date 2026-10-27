# Add a docked ligand

```{index} !ligand, !docked ligand, !add a ligand, !substrate, !small molecule
```

Add a small molecule, such as a docked substrate or a crystal ligand, to a
protein simulation, and keep it in the binding site.

The shipped example `examples/ligand_restraint` is a complete system: T4
lysozyme L99A with benzene in its cavity and a restraint that keeps it there.
Start from its config to try the steps below:

```bash
polyzymd study add-condition "Benzene" --config ~/polyzymd/examples/ligand_restraint/config.yaml --study my_study
```

## Step 1: Make the ligand's SDF

PolyzyMD places the ligand exactly where its SDF puts it. The SDF must hold:

- the pose, in the coordinates of your prepared protein PDB: a docking result
  for that same PDB, or the ligand of a crystal structure;
- every hydrogen, and the bond orders.

For a crystal ligand, {ref}`ligand-sdf-from-crystal` makes the SDF from the
PDB file and the ligand's SMILES. Prepare the protein first, as
{doc}`../tutorials/prepare_pdb_for_openff` shows.

## Step 2: Add the ligand to the config

Add a `small_molecule` entry under `solutes:`, after the protein:

```yaml
solutes:
  - type: protein
    # ...
  - type: small_molecule
    name: "benzene"
    sdf_path: "structures/benzene.sdf"
    residue_name: "BNZ"
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

`name` is the handle that restraints and analyses use. `residue_name` is the
ligand's residue name in the built system. An SDF with several conformers
places the one that `conformer_index` names (0 by default). For every key, see
{doc}`../reference/configuration`.

## Step 3: Check the config

```bash
polyzymd validate -c config.yaml
polyzymd build -c config.yaml -r 1 --dry-run
```

`validate` checks the entry and lists the ligand in its summary, as
`Substrate: benzene`. `build --dry-run` prints what the build would make,
without making it.

## Step 4: Name the ligand's atoms

In the built system, a ligand atom is named by its element, its count among
the atoms of that element in the SDF, and an `x`: `C1x`, `C2x`, ..., `H1x`.
Write `name C1x` in a selection, not `name C1`. See
{ref}`restraint-ligand-atom-names`.

## Step 5: Keep the ligand in the site

A docked ligand can leave its site during equilibration or production. To keep
it there, add a distance restraint between a ligand atom and a protein atom of
the site. {doc}`restraints` shows how; a flat-bottom restraint applies no force
while the ligand stays within the distance you set.

Distance restraints need `engine: openmm`.

## Next steps

- {doc}`restraints`: restrain the ligand in one equilibration stage only, or
  choose another restraint type.
- {doc}`../reference/configuration`: every key of the `small_molecule` entry.
