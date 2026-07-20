# Native GLYCAM exact OpenMM and GROMACS reference

This reference describes the exact GLYCAM route selected by:

```yaml
force_field:
  glycan_policy: strict_glycam
  conjugate_parameterization: native_openmm_glycam
```

## Force-field ownership

| Domain | Owner | Notes |
|--------|-------|-------|
| Protein | Amber ff14SB | Loaded through OpenMM's Amber14 XML stack. |
| Modified Asn | GLYCAM `NLN` | PolyzyMD maps each modified Asn/ASX residue to `NLN` in `residueTemplates`. |
| Glycan residues | GLYCAM_06j-1 | Input PDBs must use canonical GLYCAM names and complete `CONECT` records. |
| Water | TIP3P | Native route currently accepts TIP3P water. |
| Disconnected DMSO/free polymer/organics | OpenFF Sage through `SMIRNOFFTemplateGenerator` | Component must be disconnected from Amber/GLYCAM domains and have assigned partial charges. Disconnected multi-residue Sage polymers are represented as one OpenMM residue solely for template matching; original monomer segmentation remains audited. |

The backward-compatible default route is Sage/Interchange:

```yaml
force_field:
  glycan_policy: sage_fallback
  conjugate_parameterization: openff_interchange
```

`strict_glycam` is an explicit opt-in and fails closed. It does not silently use
OpenFF Sage for glycan residues.

## Pablo scoped identities and canonical-name restoration

PolyzyMD writes a full modified topology for Pablo parsing. Repeated or branched
GLYCAM residues can otherwise be ambiguous when a product contains repeated
canonical residue names. To keep matching deterministic, PolyzyMD assigns
attachment-scoped internal residue aliases for the Pablo-only ingestion PDB,
records the canonical residue name and scoped alias in attachment metadata, and
restores canonical GLYCAM residue names and atom names on the OpenFF topology
after Pablo succeeds.

The restored topology carries metadata such as `canonical_residue_name`,
`pablo_scoped_residue_name`, and `product_identity_source` for aliased residues.
The native OpenMM handoff then maps each modified Asn residue to GLYCAM `NLN`.
Multi-glycan systems such as Asn25 plus Asn60 use the same scoped identity model
for each attachment.

## ExactExportBundle

The native GLYCAM route returns an `ExactExportBundle`, not a public vanilla
Interchange. The bundle contains:

- authoritative OpenMM `System`, `Topology`, and positions;
- a private baseline Interchange used only by exporter internals;
- `exact_openmm_exceptions.json`, a schema v2 sidecar containing authoritative
  OpenMM atom identities, topology bonds, `NonbondedForce` particles,
  exceptions, constraints, route invariants, and hashes;
- `native_openmm_glycam_audit.json`, the route audit.

The wrapper exists because OpenFF Interchange 0.5.x models 1-4 scaling as a
global `scale14` convention. GLYCAM mixed scaled/unscaled explicit 1-4 exceptions
are local OpenMM exception records with per-exception provenance. A vanilla
Interchange export cannot represent that table losslessly.

The current limitation is upstream representation parity. PolyzyMD keeps the
native OpenMM system and versioned sidecar as the authoritative export contract
until vanilla Interchange can import and export explicit OpenMM
`NonbondedForce` exceptions without losing local exception semantics.

## Meaning of exact

For this route, **exact** means:

- force-field parameter ownership is explicit: ff14SB protein, GLYCAM glycan and
  `NLN`, TIP3P water, and disconnected precharged Sage components;
- every local OpenMM exception/exclusion needed by the GLYCAM route is captured
  in the sidecar and audited before GROMACS patching;
- the GROMACS topology patch fails closed on count, order, topology, hash, pair,
  or exclusion mismatch.

It does **not** mean bitwise-identical full PME total energies between OpenMM and
GROMACS. PME mesh spacing, interpolation order, Ewald tolerance, switching or
modifier choices, and related engine controls remain user-controlled simulation
hyperparameters. Configure those settings deliberately and validate short
engine-parity or ensemble checks for your production protocol. The audit gates
address local parameters and exception/exclusion transfer; they do not remove
normal PME implementation differences between engines.

## Native route invariants

The native OpenMM GLYCAM system is created with fixed invariants that are recorded
in the native GLYCAM audit and sidecar:

| Setting | Value |
|---------|-------|
| Nonbonded method | PME |
| Nonbonded cutoff | 1.0 nm |
| Constraints | `HBonds` |
| Rigid water | `true` |

## OpenMM execution semantics

OpenMM callers use the native objects directly:

- `ExactExportBundle.to_openmm()` returns the authoritative OpenMM `System`.
- `ExactExportBundle.to_openmm_topology()` returns the authoritative OpenMM
  `Topology`.
- `ExactExportBundle.to_openmm_positions()` returns authoritative positions.

The private baseline Interchange is not the runtime authority for OpenMM.

## Exact GROMACS export semantics

PolyzyMD's GROMACS exporter recognizes an exact bundle and uses this sequence:

1. Ask the private baseline Interchange to write baseline `.gro`/`.top` files.
2. Remove the Interchange stub MDP.
3. Preflight count, atom-order identity, topology-bond identity, and normalized
   GROMACS output hashes against `exact_openmm_exceptions.json`; fail closed
   before patching on mismatch.
4. Patch the topology using `exact_openmm_exceptions.json`.
5. Set `gen-pairs` to `no` in `[ defaults ]`.
6. Set molecule-type `nrexcl` to `0` so automatic bonded exclusions are not the
   authority.
7. Replace local `[ pairs ]` sections with explicit function-2 rows for every
   nonzero OpenMM exception.
8. Replace local `[ exclusions ]` sections with exact zero and nonzero exception
   exclusions.
9. Re-expand local rows back to global atom pairs and fail closed if the patched
   topology differs from the sidecar.
10. Write `<prefix>_exact_gromacs_audit.json` with counts, hashes, and mismatch
   counters.

Raw `Interchange.to_gromacs()` on the underlying private Interchange must not be
used for this route.

## Provenance files

| File | Purpose | Key fields |
|------|---------|------------|
| `native_openmm_glycam_audit.json` | Native GLYCAM route audit | `route`, `residue_templates`, `glycam_template_matches`, `crosslinks`, `sage_template_generator`, `route_invariants`, `unsupported_boundary_diagnostics` |
| `exact_openmm_exceptions.json` | Schema v2 exact sidecar | `schema_version`, `atoms`, `topology_bonds`, `gromacs_atoms`, `gromacs_topology_bonds`, `atom_order_hash`, `topology_hash`, `gromacs_atom_order_hash`, `gromacs_topology_hash`, `route_invariants`, `exception_hash`, `provenance` |
| `<prefix>_exact_gromacs_audit.json` | GROMACS patch validation | `route`, raw/patched pair and exclusion counts, mismatch counts, `exception_hash`, `atom_order_hash`, `gromacs_atom_order_hash`, `gromacs_topology_hash`, `route_invariants` |

Use `python -m json.tool <file>` or `jq` to inspect these files.

## Limitations

- Canonical GLYCAM-named glycan PDBs with complete `CONECT` records are required.
- Covalently attached Sage polymer across an Amber/GLYCAM boundary is unsupported;
  disconnected precharged Sage components are supported.
- Disconnected multi-residue Sage polymers may appear as one OpenMM residue in
  the native route only for `SMIRNOFFTemplateGenerator` matching. Use audit
  provenance, not that internal residue, to recover monomer identities.
- Automatic ion placement and neutralization are not audited for this route.
- Exact GROMACS export currently rejects `component_info` position-restraint
  postprocessing. Use the exact OpenMM route with staged equilibration
  `position_restraints`, or use the explicit Sage/OpenFF route when GROMACS
  position-restraint postprocessing is required.
- Pablo must be able to represent the product topology. Degree, residue-template,
  or repeated-residue matching failures are blockers, not warnings.
- Only PolyzyMD exact OpenMM and exact GROMACS exports are guaranteed. Arbitrary
  vanilla Interchange exporters are outside the guarantee.
