# Glycosylate a protein with GLYCAM06j

```{index} !glycosylation, !glycan, !glycam06j, !n-glycan, !o-glycan
```

Put N- or O-linked glycans on a protein and give them GLYCAM06j parameters.
For a first walk-through, do {doc}`../tutorials/glycosylated_protein`.

## Step 1: Get the glycan files

PolyzyMD reads each glycan from a PDB file with GLYCAM06 residue and atom
names. Build one with [GLYCAM-Web](https://glycam.org/) or download one from
[GlyGen](https://www.glygen.org/) by its GlyTouCan accession, and save it as
`structures/glycans/glycam_<accession>.pdb`. That folder is the default
`glycan_library`.

Each file needs:

- **GLYCAM06 residue names** (`0YB`, `4YB`, `VMB`, ...). GLYCAM06j finds its
  templates by them.
- **The `ROH` reducing end** as residue 1. Its O1 and HO1 leave when the
  first sugar bonds to the protein. A file without `ROH` is refused: its
  anomeric carbon has nothing to leave.
- **`CONECT` records for every bond.** GlyGen's files have none. Add them in
  PyMOL (load the file, `save` it with `set pdb_conect_all, on`) and check the
  result: PolyzyMD compares every glycan residue's bonds with its GLYCAM
  template and names any missing bond.

## Step 2: Name the sites

A site is `chain:residue`, such as `A:N60` or `A:ASN60`. Write the glycans
as a map from accession to residues:

```yaml
    modifications:
      - attach:
          type: glycan
          reaction: N-glycosylation
          chains: [A]
          glycans:
            G42666HT: [N60]
            G55220VL: [N25]
    modification_parameters:
      glycan: {force_field: [GLYCAM_06j-1.xml], charges: {method: library}}
```

- `reaction` decides which atoms bond and leave: `N-glycosylation` bonds C1
  to the Asn ND2 (HD22 and the `ROH` group leave); `O-glycosylation` bonds C1
  to a Ser OG or Thr OG1.
- The bonding atoms default, so you rarely write them: the protein atom is
  Asn ND2, Ser OG or Thr OG1 (`atom_name`), and the glycan atom is C1
  (`fragment_atom_name`) of the first residue in the file that has one
  (`fragment_resnum`; residue 1 of a GLYCAM-Web or GlyGen file is the `ROH`
  aglycone, so this is the reducing-end sugar). Give any of the three to
  override it; the build record lists the atoms used.
- `chains: [A, C, E]` puts the same glycans on every listed chain, for a
  homo-oligomer such as the BG505 SOSIP trimer.

The BG505 trimer, with 60 gp120/gp41 glycans and 3 Fab glycans, is three
entries: one for chains A, C and E, one for B, D and F, and one for the Fab
light chains.

### Sites from a table

When the sites come from a spreadsheet, list them in a CSV file with a `site`
and a `glycan` column:

```text
site,glycan
A:N60,G42666HT
A:N25,G55220VL
```

```yaml
      - attach:
          type: glycan
          reaction: N-glycosylation
          sites_file: structures/glycan_sites.csv
```

Both forms expand to the same records, one per site.

## Step 3: Check the plan

```bash
polyzymd validate -c config.yaml --plans
```

Every site gets a row: its segment (`G001`, `G002`, ...), the residue, the
glycan file and reaction, and the atoms added and removed. Add
`--format agent` for JSON.

## Step 4: Build

```bash
polyzymd build -c config.yaml
```

The physics audit checks each glycan's parameters against OpenMM's own
GLYCAM06j, including the 1-4 pairs that GLYCAM06j leaves unscaled.

## What is refused, and why

- **A glycan on GLYCAM06j on the first or last residue of a chain.** The
  GLYCAM linkage residue needs a neighbour on each side. If you meant to do this, please
  [open an issue](https://github.com/joelaforet/polyzymd/issues).
- **A residue name GLYCAM06j does not define.** Rename the residues, or
  rebuild the glycan with GLYCAM-Web.
- **A glycan on another force field** builds, with the warning
  `parameters.glycan_force_field`. GLYCAM06j is recommended: a comparison with
  Sage 2.3.0 found similar charges but stiffer glycosidic angles and different
  torsion minima (see {doc}`../explanation/modified_residues`). To use another
  force field on purpose, accept the warning in `modeling_exceptions` with
  your reason.
