# Analyze a modified protein

```{index} !segments, !segid, !modified protein analysis
```

Select a modified protein and its attachments in analyses.

## The segments

A protein with its attachments is one component, written in several
segments:

| Segment | Holds |
|---|---|
| `PROT` | The protein's residues. Mutated and modified residues stay here, as do groups merged into a residue (`merge: true`) |
| `G001`, `G002`, ... | One glycan each |
| `P001`, ... | One grafted polymer chain each |
| `M001`, ... | One other attached group each |

Residues of an attachment segment are numbered from 1, so `segid G001 and
resid 1` is the first sugar of the first glycan. The segment is in columns
73 to 76 of the PDB files, in the PSF, and in the GROMACS index groups.

## Select the protein

MDAnalysis's `protein` keyword matches canonical residue names only, so it
misses non-canonical residues (AzF, a renamed glycosylated Asn) and would
silently drop them from an RMSD or RMSF. Select the segment instead:

```python
protein = universe.select_atoms("segid PROT")
glycans = universe.select_atoms("segid G001 G002")
```

In PyMOL, use `segi PROT`.

Residue numbers repeat across segments: `resid 1` is the protein's first
residue and also the first residue of every glycan and grafted chain. Pair a
residue number with its segment, as in `segid PROT and resid 60` or
`segid G001 and resid 1`. `polyzymd validate` and `polyzymd study check` warn
when a restraint or analysis selection names residues by number alone and
would also match atoms of an attachment segment.

## Contacts with grafted and soluble chains

The `contacts` analysis takes `polymers: soluble`, `attached` or `all`.
`soluble` (the default) measures the packed chains, `attached` the grafted
chains, and `all` both. With a null `protein_selection`, the protein is the
`PROT` segment and attachments are left out.

## Run folders

The run folder of a modified protein ends its protein name with `-mod` and a
six-character (hex) hash of the declared changes, so a modified and an unmodified
protein of the same name never share a folder. A polymer with a `sequence`
block adds `-seq` and a hash the same way.
A config without these blocks keeps the folder name it had in PolyzyMD 1.3.
