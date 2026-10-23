# Rewrite a config for 1.3.0rc5

PolyzyMD 1.3.0rc5 changes the layout of `config.yaml` once, so that every
component of the system states its own force field and charges. A config
written for 1.3.0rc4 or earlier stops with an error that links here. Rewrite
it by hand, or start again from `polyzymd init` and copy your values in.
Later 1.3 releases only add keys; an rc5 config keeps working.

## What moved

| 1.3.0rc4 and earlier | 1.3.0rc5 |
|---|---|
| `enzyme:` | an entry of `solutes:` with `type: protein` |
| `substrate:` | an entry of `solutes:` with `type: small_molecule` |
| `polymers:` | an entry of `solutes:` with `type: polymer` |
| `polymers: {enabled: false}` or `polymers: null` | leave the polymer entry out |
| `force_field.protein` | `force_field:` on the protein entry, as a list |
| `force_field.small_molecule` | `force_field:` on `solvent:` and on every other solute (the same file) |
| `force_field.nagl_model` | `model:` in each `charges: {method: nagl}` |
| `substrate.charge_method: nagl` | `charges: {method: nagl}` |
| `polymers.charger` (dynamic polymers) | `charges: {method: nagl}` (or `am1bcc`, `espaloma`) |
| cached polymers (charges in the SDFs) | `charges: {method: preset}` |
| `co_solvents[].charge_method` | `co_solvents[].charges` |
| `solvent.primary: {type: water, model: tip3p}` | `solvent.water: {model: tip3p}` |
| `restraints[].atom1: {selection: "protein and resid 105 and name OG"}` | `atom1: {component: <protein name>, selection: "resid 105 and name OG"}` |

A config that did not name a force field was built with
`openff-2.0.0.offxml` for everything but the protein. To build the same
physics, write that file; the templates now show `openff-2.3.0.offxml`.

## Example

Before:

```yaml
enzyme:
  name: "LipA"
  pdb_path: "structures/LipA_prepared.pdb"

substrate:
  name: "Resorufin-Butyrate"
  sdf_path: "structures/resorufin_butyrate_docked.sdf"
  charge_method: "nagl"
  residue_name: "RBY"

solvent:
  primary:
    type: "water"
    model: "tip3p"
  ions:
    nacl_concentration: 0.1

restraints:
  - type: "flat_bottom"
    name: "substrate_active_site"
    atom1:
      selection: "protein and resid 105 and name OG"
    atom2:
      selection: "resname RBY and name C1"
    distance: 3.3

force_field:
  small_molecule: "openff-2.3.0.offxml"
```

After:

```yaml
solutes:
  - type: protein
    name: "LipA"
    pdb_path: "structures/LipA_prepared.pdb"
    force_field: ["ff14sb_off_impropers_0.0.4.offxml"]
  - type: small_molecule
    name: "Resorufin-Butyrate"
    sdf_path: "structures/resorufin_butyrate_docked.sdf"
    residue_name: "RBY"
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}

solvent:
  force_field: ["openff-2.3.0.offxml"]
  water:
    model: "tip3p"
  ions:
    nacl_concentration: 0.1

restraints:
  - type: "flat_bottom"
    name: "substrate_active_site"
    atom1:
      component: "LipA"
      selection: "resid 105 and name OG"
    atom2:
      component: "Resorufin-Butyrate"
      selection: "name C1"
    distance: 3.3
```

## New keys

- **`type`** (required) picks the component class: which keys the entry
  takes and how its physics is checked. It is never guessed.
- **`name`** (optional) is the handle restraints and analysis use. It
  defaults to the type. Keep your old enzyme and substrate names: run folder
  names (`{enzyme}`, `{substrate}`) and analysis identities are built from
  them, so the same names find your existing runs.
- **`segid`** (optional) is the segment ID of the component's atoms in the
  built PDB (columns 73-76) and in analysis: `PROT`, `LIG` and `POLY` by
  default, `SOL` for water and co-solvents and `ION` for ions. Select a
  component with `segid PROT` in MDAnalysis or `segi PROT` in PyMOL. GROMACS
  builds also write `components.ndx` with one group per component.
- **`charges`** names where partial charges come from: `library` (the
  protein force field's own), `preset` (the charges stored in the SDF files),
  `nagl` (with `model`, by default `openff-gnn-am1bcc-1.0.0.pt`), `am1bcc`
  or `espaloma`. Co-solvents PolyzyMD ships charges for, such as DMSO, use
  them unless you name another method.

## Limits of this release

- One solute of each type, and a protein is required.
- One force field file per solute; every solute but the protein, and the
  solvent, use the same file.
- One NAGL model for every component charged with `nagl`.

PolyzyMD tells you when a config goes beyond these, rather than building
something else.
