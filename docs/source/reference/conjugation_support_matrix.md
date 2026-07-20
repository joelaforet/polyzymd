# Conjugation Support Matrix

This reference summarizes what the current conjugation implementation supports,
what is validated, and what remains planned. It is intentionally conservative:
validation artifacts improve reliability, but they do not turn arbitrary reaction
chemistry into a supported workflow.

## Support levels

| Level | Meaning |
|-------|---------|
| Planned | Vocabulary or design intent exists, but users should not expect executable support. |
| Wired | Code paths or tests exercise part of the workflow, usually for mechanism extension. |
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
| SMILES moiety with `mechanism.name: n_glycosylation` | Wired, experimental | The path uses the generic moiety provider, resolved attachment plans, attachment specs, product-state charge patching, and validation hooks. It is not the same validated vertical slice as NHS-Lys polymer conjugation. |
| Residue-resolved PDB-fragment input through `moiety.input_path` | Executable loader; mechanism-gated assembly | The generic loader validates one single-model, serial-safe, connected PDB fragment with complete curated `CONECT` records and residue mapping. PolyzyMD rejects missing and detectably invalid explicit graphs, including unknown endpoints, self-bonds, disconnected or isolated atoms, invalid explicit-H degree, and obvious overvalence where applicable. It cannot prove every expected chemical bond is present and never repairs or infers omitted connectivity from coordinates. RDKit assigns bond orders only on the accepted exact connectivity. Mechanisms decide compatibility. N-glycosylation accepts residue-resolved multi-residue glycan PDB fragments that contain a structurally valid reducing-end anomeric C1 with explicit hydroxyl O/H atoms, including the supported separate-residue `ROH` hydroxyl-cap convention. |
| Strict native GLYCAM N-glycosylation with `force_field.conjugate_parameterization: native_openmm_glycam` | Validated vertical slice for the stated route | Requires explicit opt-in with `glycan_policy: strict_glycam`, canonical GLYCAM-named glycan PDBs with complete `CONECT`, `moiety.force_field_domain: glycan`, TIP3P water, and no automatic ion placement. PolyzyMD attaches to Asn, uses scoped internal Pablo aliases for repeated or branched residue matching, restores canonical names, maps modified Asn to GLYCAM `NLN`, builds an authoritative native OpenMM System from ff14SB + GLYCAM_06j-1 + TIP3P with PME, 1.0 nm cutoff, `HBonds` constraints, and rigid water, then writes exact OpenMM/GROMACS sidecars. Exact means local OpenMM exception/exclusion semantics, not bitwise full PME equality. Multi-glycan Asn25/Asn60-style systems are supported. |
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
   and a local patch charge bridge near the linkage using the pre-production
   OpenFF NAGL model `openff-gnn-am1bcc-0.1.0-rc.3.pt`. 
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
  explicit linkages are experimental paths and should be treated as
  experimental even when they produce artifacts.
- Generic PDB-fragment input is mechanism-gated. Strict native GLYCAM
  N-glycosylation is the documented production path for canonical GLYCAM-named
  glycan PDB fragments; arbitrary PDB-fragment continuations remain
  mechanism-specific.
- Existing configs that omit glycan route fields continue the backward-compatible
  Sage/Interchange route. Strict GLYCAM requires explicit `glycan_policy:
  strict_glycam` plus `conjugate_parameterization: native_openmm_glycam`.
- Raw vanilla Interchange exporters are not guaranteed for native GLYCAM exact
  bundles. Use PolyzyMD's authoritative OpenMM path or exact GROMACS exporter.
- OpenMM and GROMACS PME mesh/order/tolerance/modifier settings are engine
  hyperparameters. Configure and validate them deliberately; differences there do
  not imply a PolyzyMD local parameter or exception-transfer mismatch.
- Disconnected multi-residue Sage polymers are internally collapsed to one OpenMM
  residue only for Sage template matching, with monomer provenance retained in
  audits. Covalently attached Sage polymer across Amber/GLYCAM remains
  unsupported.
- Exact GROMACS export currently rejects `component_info` position-restraint
  postprocessing. Use exact OpenMM staged `position_restraints` or the explicit
  Sage/OpenFF route when GROMACS position restraints are required.
- Platform fail-fast behavior for CUDA/OpenMM incompatibilities is a follow-up;
  validate the pixi environment and OpenMM platform before submitting production
  jobs.
