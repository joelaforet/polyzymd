# Modified proteins

Ubiquitin (PDB 1UBQ, protonated at pH 7) with a covalent modification, in
TIP3P water with 0.15 M NaCl.

| Config | Change | Parameters of the changed atoms |
|---|---|---|
| `config_glycan.yaml` | chitobiose (GlyGen G42666HT) N-linked to Asn60 | GLYCAM06j |
| `config_octanoyl.yaml` | Lys63 deprotonated, then an octanoyl amide on NZ | Sage 2.3.0, NAGL charges |

The tutorial "From an idea and a PDB to a glycosylated protein" walks through
the first one.

## Run

In the `build` environment, from a copy of this folder:

```bash
polyzymd validate -c config_glycan.yaml --plans   # one row per site
polyzymd build -c config_glycan.yaml               # about 1 min
polyzymd run -c config_glycan.yaml                 # 200 ps
```

The configs run on an NVIDIA GPU (`openmm: platform: CUDA`), about 1 minute
on an RTX 5070 Ti laptop GPU. Set `platform: "CPU"` on a computer without one.

## Files

- `structures/ubiquitin.pdb`: ubiquitin, PDB 1UBQ, with hydrogens added at
  pH 7 (PDB data are in the public domain).
- `structures/glycans/glycam_G42666HT.pdb`: chitobiose with GLYCAM06 names,
  its `ROH` reducing end and CONECT records, from GlyGen accession G42666HT
  (GlyGen data are CC BY 4.0).
