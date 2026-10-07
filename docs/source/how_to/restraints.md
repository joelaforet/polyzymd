# Add Distance Restraints

Use this guide when you want to keep two atoms near a target distance during a
simulation, such as holding a substrate near an active site.

A distance restraint acts in every phase of the run, production included. It
has no per-stage setting. To hold a substrate in one equilibration stage only,
use a position restraint in that stage instead: see
{ref}`restrain-substrate-one-stage`.

## Step 1: choose the restraint type

PolyzyMD supports four distance restraint styles:

| Type | Best for | Behavior |
|------|----------|----------|
| `flat_bottom` | most active-site restraints | no force inside the cutoff, harmonic outside |
| `harmonic` | fixed target distances | harmonic force at all distances |
| `upper_wall` | keeping atoms from drifting apart | same practical behavior as an upper bound |
| `lower_wall` | preventing atoms from getting too close | harmonic only below the cutoff |

For most ligand-placement workflows, start with `flat_bottom`.

## Step 2: add the restraint to `config.yaml`

Example:

```yaml
restraints:
  - type: "flat_bottom"
    name: "substrate_active_site"
    atom1:
      selection: "protein and resid 76 and name OG"
      description: "Catalytic serine oxygen"
    atom2:
      selection: "resname LIG and name C1"
      description: "Substrate carbonyl carbon"
    distance: 3.5
    force_constant: 10000.0
    enabled: true
```

The units differ between the two keys:

- `distance` is in Å.
- `force_constant` is in kJ mol⁻¹ nm⁻². 4184 kJ mol⁻¹ nm⁻² is
  10 kcal mol⁻¹ Å⁻².

## Step 3: make selections specific

PolyzyMD reads selections with MDTraj. Before MDTraj reads a selection,
PolyzyMD translates `resid`, `chain` and `pdbindex` to the meanings in the
table below. The most useful selectors are:

| Keyword | Meaning | Example |
|---------|---------|---------|
| `resid` | residue number | `resid 76` |
| `resname` | residue name | `resname LIG` |
| `name` | atom name | `name OG` |
| `pdbindex` | Position in the built system, counted from 1: the PDB ATOM serial PolyzyMD writes. Analyses read it the same way. | `pdbindex 2740` |
| `index` | OpenMM atom index, 0-indexed | `index 2739` |
| `chain` | chain identifier | `chain A` |

Combine them with `and`, `or`, `not` and parentheses. Write a range as
`resid 70 to 80`.

- Atom and residue names are case-sensitive: `name OG` and `name og` are
  different selections.
- The words `and`, `or`, `not` and `to` may be written in any case.
- `index`, `resid` and `pdbindex` take numbers only. A word after one of them,
  such as `index 6 x`, is refused with an error that names the selection.
- Each restraint selection must match exactly one atom. Otherwise the build
  stops with `Restraint '<name>' requires exactly one atom per selection. Got
  <n> for atom1, <m> for atom2`. A selection that matches no atom stops with
  `No atoms match selection: '<selection>'`.

```{warning}
Always make protein selections chain-aware enough to avoid accidental matches.
`protein and resid 76 and name OG` is safer than `resid 76 and name OG`.
```

```{important}
`resid` refers to the residue number in the **built** topology, not in your
input PDB. PolyzyMD renumbers the protein consecutively from 1 (chain A), so
if your PDB starts at residue 5, crystal-structure residue 144 becomes
`resid 140`. Substrate (chain B) residues also restart at 1. Check
`solvated_system.pdb` after `polyzymd build` and adjust selections and
analysis definitions accordingly.
```

## Step 4: find the right atom indices

If the atom names in your input files are not enough, first build the system so
you can inspect `solvated_system.pdb`:

```bash
pixi run -e build polyzymd build -c config.yaml
```

Then open the solvated structure in PyMOL and inspect the atom serial shown in
the built system, not only in the original ligand or protein input file.

That built PDB is the best source for `pdbindex` values because it reflects the
final atom ordering used by the simulation.

## Step 5: validate and test

Run:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd build -c config.yaml --dry-run
```

During a real build or run, PolyzyMD should report that the restraint was
applied.

## Common patterns

### Keep a substrate near the catalytic residue

```yaml
restraints:
  - type: "flat_bottom"
    name: "substrate_catalytic"
    atom1:
      selection: "protein and resid 76 and name OG"
    atom2:
      selection: "resname LIG and name C1"
    distance: 3.5
    force_constant: 10000.0
    enabled: true
```

(restrain-substrate-one-stage)=
### Hold a substrate during the first equilibration stage only

Position restraints are set per equilibration stage. The `ligand_heavy` group
holds the substrate heavy atoms at their starting coordinates. List it in the
first stage only, and leave it out of the later stages:

```yaml
simulation_phases:
  equilibration_stages:
    - name: "restrained"
      # ... duration, ensemble and temperature of the stage
      position_restraints:
        - group: "protein_heavy"
          force_constant: 4184.0
        - group: "ligand_heavy"
          force_constant: 4184.0
    - name: "free"
      # ... no position_restraints: the substrate moves freely
```

Production has no position restraints. The OpenMM run log prints
`Removing N position restraint force(s) for next stage` when a stage ends.
See {doc}`equilibration` for the other groups.

### Restrain a protein-protein distance

```yaml
restraints:
  - type: "harmonic"
    name: "domain_distance"
    atom1:
      selection: "protein and resid 50 and name CA"
    atom2:
      selection: "protein and resid 150 and name CA"
    distance: 25.0
    force_constant: 100.0
    enabled: true
```

### Temporarily disable a restraint

```yaml
restraints:
  - type: "flat_bottom"
    name: "optional_restraint"
    atom1:
      selection: "protein and resid 76 and name OG"
    atom2:
      selection: "resname LIG and name C1"
    distance: 4.0
    force_constant: 5000.0
    enabled: false
```

## Force constant starting points

| Use case | Suggested `force_constant` (kJ mol⁻¹ nm⁻²) |
|----------|----------------------------|
| strong restraint | `10000-50000` |
| moderate restraint | `1000-5000` |
| weak guiding restraint | `100-500` |

## Troubleshooting

### no atoms match the selection

Check residue numbering, atom names, and chain identity in the built PDB.
Remember that the protein is renumbered from 1: a PDB whose first residue is
number 5 shifts every `resid` down by 4 relative to the crystal numbering.

### selection matches more than one atom

Make the selection more specific by adding `name`, `chain`, or `pdbindex`.

### restraint seems to do nothing

Check that:

- `enabled: true` is set
- the distance is in angstroms
- the force constant is large enough for your use case

## Related pages

- schema details: {doc}`../reference/configuration`
- staged setup workflows: {doc}`equilibration`
- first build tutorial: {doc}`../get_started/quickstart`

<!-- IMAGE OPPORTUNITY: Add a PyMOL screenshot of `solvated_system.pdb` with one
protein atom and one ligand atom labeled, plus an annotation showing how those
labels map to `selection` fields in YAML. -->
