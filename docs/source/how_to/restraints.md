# Add Distance Restraints

Use this guide when you want to keep two atoms near a target distance during a
simulation, such as holding a substrate near an active site during early
equilibration or production.

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
      selection: "protein and resid 77 and name OG"
      description: "Catalytic serine oxygen"
    atom2:
      selection: "resname LIG and name C1"
      description: "Substrate carbonyl carbon"
    distance: 3.5
    force_constant: 10000.0
    enabled: true
```

## Step 3: make selections specific

PolyzyMD uses MDAnalysis-style selections. The most useful selectors are:

| Keyword | Meaning | Example |
|---------|---------|---------|
| `resid` | residue number | `resid 77` |
| `resname` | residue name | `resname LIG` |
| `name` | atom name | `name OG` |
| `pdbindex` | PDB atom serial, 1-indexed | `pdbindex 2740` |
| `index` | OpenMM atom index, 0-indexed | `index 2739` |
| `chain` | chain identifier | `chain A` |

Combine them with `and` or `or` as needed.

```{warning}
Always make protein selections chain-aware enough to avoid accidental matches.
`protein and resid 77 and name OG` is safer than `resid 77 and name OG`.
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
      selection: "protein and resid 77 and name OG"
    atom2:
      selection: "resname LIG and name C1"
    distance: 3.5
    force_constant: 10000.0
    enabled: true
```

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
      selection: "protein and resid 77 and name OG"
    atom2:
      selection: "resname LIG and name C1"
    distance: 4.0
    force_constant: 5000.0
    enabled: false
```

## Force constant starting points

| Use case | Suggested `force_constant` |
|----------|----------------------------|
| strong restraint | `10000-50000` |
| moderate restraint | `1000-5000` |
| weak guiding restraint | `100-500` |

## Troubleshooting

### no atoms match the selection

Check residue numbering, atom names, and chain identity in the built PDB.

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
- first build tutorial: {doc}`../tutorials/quickstart`

<!-- IMAGE OPPORTUNITY: Add a PyMOL screenshot of `solvated_system.pdb` with one
protein atom and one ligand atom labeled, plus an annotation showing how those
labels map to `selection` fields in YAML. -->
