# Build an N-glycosylated coordinate artifact from a GlyGen PDB

Use this guide when an RCSB PDB entry identifies the glycan you want to model,
and a GlyGen or GlyCAM-style source supplies a downloadable residue-resolved 3D
PDB that PolyzyMD can attach to an asparagine site with the current PDB-fragment
N-glycosylation workflow.

```{important}
This workflow is **coordinate-only by default**. It writes a residue-resolved
crosslinked PDB plus provenance, but it does not perform GLYCAM, CHARMM, OpenFF,
or Pablo parameterization. Treat the output as an inspected structural handoff
for external workflows.
```

## Prerequisites

- A cleaned protein PDB containing the target Asn with explicit hydrogens.
- A GlyGen or GlyCAM-style glycan PDB for the desired glycan. GlycoShape PDBs
  may work only after reducing them to one model and confirming that they contain
  the strict `ROH` residue with `O1`/`HO1` leaving-atom markers described below.
- A PolyzyMD configuration that enables `conjugation` and defines one
  `n_glycosylation` attachment using `moiety.input_path`.

The current PDB-fragment ingestion MVP supports exactly one GlyGen/GlyCAM-profile
glycan attachment per coordinate-only build and does not mix PDB-fragment inputs
with SMILES or polymer moiety sources.

## 1. Find the glycan on RCSB PDB

1. Open the RCSB PDB structure page, for example
   <https://www.rcsb.org/structure/5FYJ>.
2. Scroll to the **Glycans** section.
3. Identify the glycan you want to model and note its GlyGen accession, such as
   `G80966KZ`.

## 2. Download a 3D glycan PDB from a GlyGen/GlyCAM-style source

1. Open the GlyGen glycan page, for example
   <https://www.glygen.org/glycan/G80966KZ>.
2. Use the GlyGen accession to find a downloadable 3D PDB from GlyGen itself or
   a linked GlyCAM-style source. RCSB identifies the glycan in the protein
   structure; it does not by itself provide the standalone reducing-end glycan
   PDB that this loader expects.
3. GlycoShape downloads may be usable only after manual preparation: extract a
   single `MODEL` and verify the strict `ROH` residue with atoms named `O1` and
   `HO1` is present.
4. Save it in your project, for example:

   ```text
   structures/G80966KZ_glycam.pdb
   ```

### File considerations

- The generic loader requires a single connected fragment graph and preserves
  residue labels.
- Files with `CONECT` records use those records directly.
- Files without `CONECT` records use RDKit coordinate-inferred connectivity and
  record `coordinate_inferred` provenance in the sidecar. Inspect the output PDB
  and sidecar before using the structure downstream.
- Multi-model PDB files are not automatically reduced. If a GlycoShape or other
  download contains more than one `MODEL` record, extract the model you want into
  a single-model PDB before loading it.
- The input must contain one `ROH` residue with atoms named `O1` and `HO1`. This
  marks the reducing-end leaving group used by the current N-glycosylation
  template.

## 3. Configure the N-glycosylation attachment

Add one enabled attachment under `conjugation.attachments`. Use the cleaned
protein's Asn chain and residue number, and point `moiety.input_path` at the
downloaded glycan PDB.

```yaml
enzyme:
  name: 5fyj_example
  pdb_path: structures/5fyj_clean_asn.pdb

conjugation:
  enabled: true
  attachments:
    - name: asn60_glygen_glycan
      site:
        chain_id: A
        residue_name: ASN
        residue_number: 60
        atom_name: ND2
      moiety:
        name: G80966KZ
        input_path: structures/G80966KZ_glycam.pdb
      mechanism:
        name: n_glycosylation
```

The built-in mechanism forms an Asn `ND2`--glycan reducing-end `C1` bond. It
removes one Asn `ND2` hydrogen and the glycan `ROH:O1`/`ROH:HO1` leaving atoms,
renames the protein-site residue to `ASX`, and preserves the glycan residue
labels in chain `C` after removing the `ROH` leaving group.

## 4. Run the coordinate-only build through the Python API

Use the public config-driven conjugation API. The default
`pdb_fragment_output_mode` is `coordinate_only`; the explicit settings below make
that choice visible.

```python
from pathlib import Path

from polyzymd.builders.conjugation import build_conjugate_from_config
from polyzymd.builders.conjugation import ConjugatedPolymerSystemSettings

result = build_conjugate_from_config(
    "config.yaml",
    output_dir=Path("artifacts/glygen_asn60"),
    settings=ConjugatedPolymerSystemSettings(
        pdb_fragment_output_mode="coordinate_only",
    ),
)

print(result.status)
print(result.artifact_paths["pdb_fragment_coordinate_only_pdb"])
print(result.artifact_paths["pdb_fragment_pdb_fragment_ingestion"])
```

Do not use the coordinate-only result as a ready-to-run OpenMM system. It has no
final Interchange, no `system.xml`, and no solvated topology.

## 5. Inspect the artifacts

The default coordinate-only paths are:

| Artifact | Typical path | What to check |
|----------|--------------|---------------|
| Crosslinked PDB | `artifacts/glygen_asn60/conjugate-construction/pdb_fragment_coordinate_only_conjugate.pdb` | Asn site is `ASX`; glycan residues are in chain `C`; `ROH` atoms are absent; the PDB contains linkage `CONECT`/`LINK` evidence. |
| PDB-fragment ingestion sidecar | `artifacts/glygen_asn60/conjugate-polymerist-cache/asn60_glygen_glycan_pdb_fragment_ingestion.json` | `connectivity_provenance`, residue mapping, and the `n_glycosylation_profile` with reducing `C1`, removed `ROH` atoms, and linkage diagnostics. |
| Workflow JSON | `artifacts/glygen_asn60/conjugated_polymer_system_workflow.json` | `status: coordinate_only` and artifact path metadata. |

For raw GlyGen/GlyCAM PDBs without `CONECT`, pay particular attention to the sidecar's
`connectivity_provenance: coordinate_inferred` value and linkage diagnostics.
Coordinate-inferred bonds are validated for a connected graph and plausible
inter-residue C--O glycosidic linkages, but you should still inspect the emitted
PDB and provenance before parameterizing elsewhere.

## Experimental Pablo/OpenFF continuation

`ConjugatedPolymerSystemSettings(pdb_fragment_output_mode="experimental_pablo")`
continues past the coordinate artifact into the current Pablo/OpenFF path. This
mode is experimental for GlyGen/GlyCAM-profile PDB fragments and is not a GLYCAM or CHARMM
parameterization workflow. Prefer `coordinate_only` unless you are explicitly
testing the Pablo/OpenFF continuation and are prepared to validate failures and
charge/parameter provenance yourself.

## Related reference

- {doc}`../reference/protein_modification_config`
- {doc}`../reference/conjugation_support_matrix`
