# Rewrite a 1.2 config for 1.3

PolyzyMD 1.3.0 changes the layout of `config.yaml` once, so that every
component of the system states its own force field and charges. A config
written for PolyzyMD 1.2 (such as 1.2.1), or for a 1.3.0 release candidate up
to rc4, stops with an error that links here. Rewrite
it by hand, or write a new template with
`polyzymd study add-condition LABEL --new` and copy your values in.
Later 1.3 releases only add keys; a 1.3.0 config keeps working.

## What moved

| 1.2 (and 1.3.0rc4 and earlier) | 1.3.0 |
|---|---|
| `enzyme:` | an entry of `solutes:` with `type: protein` |
| `substrate:` | an entry of `solutes:` with `type: small_molecule` |
| `polymers:` | an entry of `solutes:` with `type: polymer` |
| `polymers: {enabled: false}` or `polymers: null` | leave the polymer entry out |
| `force_field.protein` | `force_field:` on the protein entry, as a list |
| `force_field.small_molecule` | `force_field:` on `solvent:` and on every other solute (the same file) |
| `force_field.nagl_model` (or nothing: 1.2 and 1.3.0rc4 and earlier charged with `openff-gnn-am1bcc-0.1.0-rc.3.pt`) | `model:` in each `charges: {method: nagl}` |
| `substrate.charge_method: nagl` | `charges: {method: nagl, model: ...}` (see below for which model) |
| `polymers.charger` (dynamic polymers) | `charges: {method: nagl, model: ...}` (or `am1bcc`, `espaloma`) |
| `polymers.reactions` and monomer SMILES with their C=C double bond (dynamic polymers) | each monomer's repeat unit with `*` on its two linking atoms (`C=C(C)C(=O)OC` becomes `[*]CC([*])(C)C(=O)OC`), and `end_groups: {head: [...], tail: [...]}`, the caps each end can get with one `*` each and their probabilities; see {doc}`dynamic_polymers` |
| cached polymers (charges in the SDFs) | `charges: {method: preset}` |
| `co_solvents[].charge_method`, or nothing, on a co-solvent PolyzyMD ships charges for (dmso, urea, ...) | `charges: {method: preset}`: those always got their shipped charges |
| `co_solvents[].charge_method` on any other co-solvent | `co_solvents[].charges: {method: nagl, model: ...}` (or `am1bcc`, `espaloma`) |
| `solvent.primary: {type: water, model: tip3p}`, or nothing (the default was tip3p) | `solvent.water: {model: tip3p}` |
| the protein's charges (always ff14SB's own) | `charges: {method: library}` on the protein entry |
| `restraints[].atom1: {selection: "protein and resid 105 and name OG"}` | `atom1: {component: <protein name>, selection: "resid 105 and name OG"}` |

A config that did not name them was built with `openff-2.0.0.offxml` for
everything but the protein, and charged every `nagl` molecule with
`openff-gnn-am1bcc-0.1.0-rc.3.pt`. To build the same physics, write those
two names. The templates now show `openff-2.3.0.offxml` and
`openff-gnn-am1bcc-1.0.0.pt` (the model Sage 2.3.0 names); moving to them is
a change of force field and charges, so make it on purpose.

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
    charges: {method: library}
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
      selection: "name C13x"  # ester carbonyl carbon of the docked RBY
    distance: 3.3
```

The atom name changes too. This SDF carries no atom names, so the built atoms
are named by element, count in SDF order and `x`. The ester carbonyl carbon is
`C13x`. See {ref}`restraint-ligand-atom-names`.

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
- **`charges`** (required on every solute and co-solvent) names where partial
  charges come from: `library` (the protein force field's own), `preset` (the
  charges stored in the SDF files, or shipped with PolyzyMD for co-solvents
  such as DMSO), `nagl` (with a required `model`, such as
  `openff-gnn-am1bcc-1.0.0.pt`, the model Sage 2.3.0 names), `am1bcc` or
  `espaloma`.
- **`force_field`** (required on every solute and on `solvent`) and
  **`solvent.water.model`** (required): PolyzyMD assumes no force field,
  charge model or water model. Each is your choice, on record in the config
  that built the system.

## Limits of this release

- One solute of each type, and a protein is required.
- One force field file per solute; every solute but the protein, and the
  solvent, use the same file.
- One NAGL model for every component charged with `nagl`.

PolyzyMD tells you when a config goes beyond these, rather than building
something else.

## Existing builds and runs

- **A replicate that has started** (it has minimization, equilibration or
  production output) continues from its own files. `submit` and the job
  scripts do not hold it to the new checks.
- **A replicate that has not started**, built by 1.2 or by 1.3.0rc4 or
  earlier, must be rebuilt. Its build records no physics audit, and `run`
  and `submit` refuse files that no passing audit vouches for. The error
  gives the command. Rebuild after you rewrite the config:
  `polyzymd build -c config.yaml -r 1`.
- **Analysis** of existing runs works as before. If you keep the old enzyme
  and substrate names as `name:`, the analysis finds the existing results.

## What the build now checks

Every build ends with the physics audit: the built system is compared, term
by term, with what the config names. The build writes
`build_provenance.json` and `build_provenance.md` (they replace
`build_manifest.json`). Open `build_provenance.md` after your first build: it
lists each component with its force field and charges, and every check with
its result. {doc}`../explanation/physics_audit` explains the checks, and
{doc}`../reference/physics_checks` lists them.

On GROMACS, the audit also reads every `.mdp` file the export wrote and
compares it with the config. Each GROMACS job script, local or SLURM, checks
`gromacs/audited_files.sha256`
before `grompp`. If you edit an exported `.top`, `.itp` or `.mdp` file by
hand, the job stops: change the config and export again instead.

## If a build fails its audit

The build stops and prints each failing check:

```text
The physics audit failed: the system is not what the config asks for.
  FAIL parameters.bonded (LipA): <which terms differ, with examples>
Details: <run folder>/build_provenance.failed.md
fix: correct the config or the inputs and build again; if a failure is a deliberate choice, accept it under modeling_exceptions with the reason (<link to the check catalog>).
```

The full record is in `build_provenance.failed.md` (for GROMACS, in
`build_provenance.md`, marked failed). Read the detail of each
failing check: it names the terms or atoms that differ from the config. Most
failures are a config that does not say what you meant, such as a force field
file you did not intend or a restraint selection that matches another atom.

A check can also warn: the build goes on and `build_provenance.md` shows
the warning. If a failure or a warning is deliberate, accept it in the config,
with the reason. For example, crystal waters kept inside a protein warn under
`geometry.buried_water`:

```yaml
modeling_exceptions:
  - check: geometry.buried_water
    component: LipA
    reason: "Crystal waters in the active-site pocket are kept on purpose."
```

The build record keeps the accepted finding and its reason. Please report any
failure you think is wrong, with the `build_provenance.failed.md` file, in an
issue on GitHub.
