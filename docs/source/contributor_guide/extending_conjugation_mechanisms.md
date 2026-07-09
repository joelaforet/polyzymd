# Extend Conjugation Mechanisms

Use this contributor guide when you are adding or reviewing a new covalent
conjugation mechanism. The current internals are still evolving, so this page is
a checklist and design contract, not a promise of a stable public plugin API.

For user-facing support levels, update {doc}`../reference/conjugation_support_matrix`.
For architecture rationale, read {doc}`../explanation/conjugation_architecture`.

## Start from a concrete mechanism

Before editing code, define the smallest realistic product you need to support.
Do not start by making arbitrary SMARTS configurable. A mechanism must specify
enough information for construction, parameterization, charge auditing, and
validation.

At minimum, write down:

- the input protein residue and atom that react;
- the moiety residue or polymer monomer that reacts;
- the product bond and expected bond order or target distance;
- all atoms removed during product formation;
- product residue names for the modified site and moiety;
- product-state atom/residue mappings after leaving atoms are removed;
- force-field or template expectations;
- charge patch and reconciliation expectations;
- validation evidence the workflow should write.

## Mechanism checklist

New mechanisms must define or derive the following pieces.

| Requirement | What to specify |
|-------------|-----------------|
| Reactive atoms | The protein-side link atom and moiety-side link atom, including residue and atom-name expectations. |
| Leaving atoms | Protein-side and moiety-side atoms that must be absent in the product. |
| Product residue names | Site and moiety residue names after reaction, including any built-in aliases. |
| Product mappings | How source atoms map to product atoms after residue renaming, atom deletion, or generated-fragment changes. |
| Product-state residue definitions | Pablo/CCD residue definitions or generated product-state libraries needed for ingestion. |
| Charge patch expectations | Which atoms keep source/template charges, which local neighborhood is patched, and what reconciliation tolerance is expected. |
| Parameter expectations | Which force fields or templates cover the final product and how missing coverage should fail. |
| Geometry expectations | Expected linkage distance range and any known relaxation or validation needs. |
| Evidence paths | JSON sidecars the validation report can audit. |

## Validation checklist

Every new mechanism or broadened support level should include tests and docs for
the canonical validation report:

- `bond_graph` passes when the expected product linkage is present and fails when
  it is missing.
- `atom_presence` detects missing link atoms and lingering leaving atoms.
- `valence_sanity` gives conservative warnings for implausible linkage bond
  counts.
- `charge_audit` records the charge bridge path, local reconciliation path,
  total/formal charge, and normalization corrections.
- `parameter_coverage` records production OpenMM system particle-count evidence
  when available.
- `linkage_geometry` records linkage distances and close-contact warnings.
- `relaxation_evidence` records `conjugate_relaxation.json` success, finite
  energies, restrained protein displacement, and linkage distance errors.

Tests should include at least one passing product and one targeted failing case
for each mechanism-specific identity or leaving-atom rule. Prefer small unit
fixtures for mechanism metadata, then one integration-style vertical slice for
the supported workflow.

## Documentation checklist

When changing support, update the docs in the same change set:

1. Update {doc}`../reference/conjugation_support_matrix` with the mechanism's
   support level and caveats.
2. Update {doc}`../reference/protein_modification_config` if config fields or
   examples change.
3. Update {doc}`../how_to/validate_conjugates` if the validation report gains or
   changes evidence fields.
4. Update tutorials only when the workflow is appropriate for first-time users.

## What not to promise

Avoid user-facing claims that exceed the implementation:

- Do not describe arbitrary chemistry as supported because a SMARTS pattern can
  be written.
- Do not claim production readiness without product-state residue definitions,
  charge evidence, parameter coverage, validation tests, and scientific caveats.
- Do not promise a stable public mechanism-plugin API while internals are still
  changing.
- Do not imply that the CUDA 12.4 simulation-only handoff is automatic unless a
  tested command or documented workflow exists.
