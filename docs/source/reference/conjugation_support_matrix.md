# Conjugation Support Matrix

This reference summarizes what the current conjugation implementation supports,
what is validated, and what remains planned. It is intentionally conservative:
validation artifacts improve reliability, but they do not turn arbitrary reaction
chemistry into a supported workflow.

## Support levels

| Level | Meaning |
|-------|---------|
| Planned | Vocabulary or design intent exists, but users should not expect executable support. |
| Wired | Code paths or tests exercise part of the workflow, usually for mechanism development. |
| Executable | PolyzyMD can run the workflow and write construction artifacts. |
| Validated vertical slice | The workflow has validation/audit artifacts and tests for the stated scope. |
| Production-ready | Scientific guidance, templates, charges, platform handoff, and validation are mature enough for routine production use. |

## Current mechanism support

| Mechanism or feature | Current support | Notes |
|----------------------|-----------------|-------|
| Opt-in config-driven conjugation with `conjugation.enabled: true` | Executable | Runs through the regular `polyzymd build -c config.yaml` path when attachments are enabled. |
| Direct public requests with `protein_pdb_path` plus `attachments` | Executable | `ConjugationEngine` and `build_conjugate()` accept the same attachment specification model used by config workflows. Raw molecule/topology keyword inputs remain unsupported. |
| NHS-lysine polymer attachment with `mechanism.name: nhs_lys_amide` | Validated vertical slice | Current reliability milestone target. Supports generic provider/resolved-plan/attachment-spec orchestration, Pablo/OpenFF Interchange, ff14SB plus polymer templates, local NAGL patch charge bridge, charge reconciliation, conjugate relaxation, validation report, and final solvation. |
| Multiple NHS-lysine polymer attachments | Executable | Multi-site attachment metadata is handled by the validation report and shared conjugate relaxation. Treat scientific validation as system-specific. |
| mBuild -> OpenFF polymer adapter for ACB parity slice | Validated adapter slice | `from_mbuild()` converts atomistic mBuild compounds to OpenFF molecules, and OpenFF/SDF adapters feed existing `GeneratedPolymerFragment` conjugation paths. ACB parity is proven against Polymerist for formula, graph, OpenFF assignments, NAGL charges, and OpenMM energy. Bundled-default dynamic generation now uses native linear methacrylate generation; custom `.rxn` workflows retain deprecated Polymerist compatibility behavior. |
| SMILES moiety with `mechanism.name: n_glycosylation` | Wired, experimental | The path uses the generic moiety provider, resolved attachment plans, attachment specs, product-state charge patching, and validation hooks. It is not the same validated vertical slice as NHS-Lys polymer conjugation. |
| Mixed mechanisms in one config | Wired, experimental | Shared assembly and validation can handle multiple resolved plans, but mixed-chemistry scientific validation remains mechanism-specific. |
| O-glycosylation | Planned | Vocabulary/design area only. |
| Arbitrary SMARTS-defined chemistry | Planned | SMARTS can describe roles for future mechanisms, but does not by itself provide placement, product residues, charge patches, templates, or validation. |
| Explicit PDB linkage | Wired | Requires explicit atom-level product residue and leaving-atom information. Treat as advanced/exploratory and validate carefully. |

## Validated NHS-Lys polymer vertical slice

The validated current path is:

1. Start from an unmodified protein PDB.
2. Enable conjugation in config and request one or more NHS-Lys polymer
   attachments.
3. Resolve the moiety through the generic provider layer from a polymer recipe,
   SMILES moiety, or explicit source supported by the mechanism.
4. Resolve the attachment plan, reactive atoms, leaving atoms, product residue
   names, and product-state mappings into an attachment spec.
5. Build a crosslinked product PDB with expected linkage `CONECT` records.
6. Ingest the product state with Pablo and build the final OpenFF Interchange.
7. Combine ff14SB-style protein treatment, polymer template charges/parameters,
   and a local NAGL/AshGC patch charge bridge near the linkage.
8. Reconcile local charge evidence instead of falling back to whole-conjugate
   AM1-BCC or unmarked cached charges.
9. Run conjugate relaxation and collect restrained OpenMM relaxation evidence when
   enabled by the build workflow.
10. Write final solvated artifacts and `conjugate_validation_report.json`.

## Validation report schema overview

The canonical report is `conjugate_validation_report.json`. It contains a
top-level `status` and these sections:

| Section | Purpose |
|---------|---------|
| `bond_graph` | Product linkage bonds expected from assembly metadata are present in product PDB `CONECT` records. |
| `atom_presence` | Link atoms are present and mechanism-declared leaving atoms are absent. |
| `valence_sanity` | Linkage atom bond counts are within conservative limits. |
| `charge_audit` | Charge bridge and local reconciliation JSON evidence are present and internally consistent. |
| `parameter_coverage` | Production OpenMM system evidence has the expected particle count when available. |
| `linkage_geometry` | Linkage distances and severe close contacts are checked in product coordinates. |
| `relaxation_evidence` | `conjugate_relaxation.json` reports successful restrained OpenMM relaxation with finite energies, stable restrained protein atoms, and checked linkage distances. |

Each section stores `checks` with a check name, status, message, and structured
evidence. See {doc}`../how_to/validate_conjugates` for a task-oriented reading
workflow.

## Artifact names

Common conjugation reliability artifacts include:

| Artifact | Meaning |
|----------|---------|
| `conjugate_validation_report.json` | Canonical validation report for the constructed product. |
| `product_state_charge_bridge.json` | Charge bridge provenance, charge totals, formal charge, and correction summary. |
| `product_state_charge_bridge_local_reconciliation.json` | Local reconciliation evidence for the patch charge bridge. |
| `conjugate_relaxation.json` | Restrained OpenMM relaxation evidence for the conjugate product. |

## Limits and caveats

- A passing validation report means available checks passed. It does not prove
  that a mechanism outside the support matrix is chemically valid.
- NHS-Lys is the validated vertical slice. Non-NHS mechanisms and generic
  explicit linkages are active development paths and should be treated as
  experimental even when they produce artifacts.
- The exact runtime handoff from a build-time OpenFF Interchange to a lean CUDA
  12.4 simulation environment is planned or workflow-specific. Follow only
  documented commands and site-tested scripts.
- Platform fail-fast behavior for CUDA/OpenMM incompatibilities is a follow-up;
  validate the pixi environment and OpenMM platform before submitting production
  jobs.
