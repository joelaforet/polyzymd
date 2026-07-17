# Build an N-glycosylated coordinate artifact from a glycan PDB

Use this guide when you have a residue-resolved multi-residue glycan PDB that
PolyzyMD can attach to an asparagine site with the current PDB-fragment
N-glycosylation workflow. GlyGen/RCSB is shown as one acquisition workflow, not
as a provenance requirement.

```{important}
This workflow is **coordinate-only by default**. It writes a residue-resolved
crosslinked PDB plus provenance, but it does not perform GLYCAM, CHARMM, OpenFF,
or Pablo parameterization. Treat the output as an inspected structural handoff
for external workflows.
```

## Prerequisites

- A cleaned protein PDB containing the target Asn with explicit hydrogens.
- A residue-resolved glycan PDB for the desired glycan. The reducing-end
  anomeric `C1` must have explicit hydroxyl O/H atoms, either as an ordinary
  residue-local hydroxyl or as the supported separate-residue `ROH` cap with
  `O1`/`HO1` atom names.
- A PolyzyMD configuration that enables `conjugation` and defines one
  `n_glycosylation` attachment using `moiety.input_path`.

The current PDB-fragment ingestion MVP supports one residue-resolved glycan PDB
attachment per coordinate-only build and does not mix PDB-fragment inputs with
SMILES or polymer moiety sources.

## 1. Find the glycan on RCSB PDB

1. Open the RCSB PDB structure page, for example
   <https://www.rcsb.org/structure/5FYJ>.
2. Scroll to the **Glycans** section.
3. Identify the glycan you want to model and note its GlyGen accession, such as
   `G80966KZ`.

## 2. Download or prepare a 3D residue-resolved glycan PDB

1. Open the GlyGen glycan page, for example
   <https://www.glygen.org/glycan/G80966KZ>.
2. Use the GlyGen accession to find a downloadable 3D PDB from GlyGen or another
   residue-resolved glycan source. RCSB identifies the glycan in the protein
   structure; it does not by itself provide the standalone reducing-end glycan
   PDB that this loader expects.
3. GlycoShape or other downloads may be usable after manual preparation: extract
   a single `MODEL`, retain residue labels, and verify explicit reducing-end
   C1-O-H chemistry.
4. Save it in your project, for example:

   ```text
   structures/G80966KZ_glycam.pdb
   ```

### File considerations

- The generic loader requires a single connected fragment graph and preserves
  residue labels.
- Files must include complete, curated `CONECT` records for the whole fragment.
  PolyzyMD rejects missing and detectably invalid explicit graphs, including
  unknown endpoints, self-bonds, disconnected or isolated atoms, invalid
  explicit-H degree, and obvious overvalence where applicable.
- PolyzyMD cannot prove every expected chemical bond is present and never
  repairs or infers omitted connectivity from coordinates. RDKit assigns bond
  orders only on the accepted exact connectivity. Inspect the output PDB and
  sidecar before using the structure downstream.
- Multi-model PDB files are not automatically reduced. If a GlycoShape or other
  download contains more than one `MODEL` record, extract the model you want into
  a single-model PDB before loading it.
- The input must contain a structurally valid reducing-end `C1` bonded to an
  explicit hydroxyl O/H group. The `ROH:O1`/`ROH:HO1` cap is accepted, but not
  required when an ordinary residue-local hydroxyl is present.

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
    - name: asn60_residue_resolved_glycan
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
removes one Asn `ND2` hydrogen and the validated glycan hydroxyl O/H leaving
atoms, renames the protein-site residue to `ASX`, and preserves the remaining
glycan residue labels in chain `C`.

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
    output_dir=Path("artifacts/nglycan_asn60"),
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
| Crosslinked PDB | `artifacts/nglycan_asn60/conjugate-construction/pdb_fragment_coordinate_only_conjugate.pdb` | Asn site is `ASX`; glycan residues are in chain `C`; hydroxyl leaving atoms are absent; the PDB contains linkage `CONECT`/`LINK` evidence. |
| PDB-fragment ingestion sidecar | `artifacts/nglycan_asn60/conjugate-polymerist-cache/asn60_residue_resolved_glycan_pdb_fragment_ingestion.json` | `connectivity_provenance`, residue mapping, and the `n_glycosylation_profile` with reducing `C1`, removed hydroxyl atoms, and linkage diagnostics. |
| Workflow JSON | `artifacts/nglycan_asn60/conjugated_polymer_system_workflow.json` | `status: coordinate_only` and artifact path metadata. |

Raw glycan PDBs must include complete curated `CONECT` records. PolyzyMD
rejects missing or detectably invalid explicit graphs, but cannot prove every
expected chemical bond is present; for accepted fragments, confirm the sidecar
reports explicit `CONECT` provenance and inspect the linkage diagnostics before
parameterizing elsewhere.

## Experimental Pablo/OpenFF continuation

`ConjugatedPolymerSystemSettings(pdb_fragment_output_mode="experimental_pablo")`
continues past the coordinate artifact into the current Pablo/OpenFF path. This
mode is experimental for residue-resolved glycan PDB fragments and is not a GLYCAM or CHARMM
parameterization workflow. Prefer `coordinate_only` unless you are explicitly
testing the Pablo/OpenFF continuation and are prepared to validate failures and
charge/parameter provenance yourself.

## Related reference

- {doc}`../reference/protein_modification_config`
- {doc}`../reference/conjugation_support_matrix`
