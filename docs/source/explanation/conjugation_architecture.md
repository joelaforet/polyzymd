# Conjugation Architecture and Validation Boundaries

This page explains the design decision behind PolyzyMD's current conjugation
workflow and the validation boundary introduced for the reliability milestone.
It is written in an ADR style so contributors can understand why the system is
conservative about supported mechanisms.

## Status

Accepted for the current conjugation reliability milestone.

## Context

PolyzyMD needs to build enzyme-polymer conjugates that can pass through protein
preparation, product-state residue handling, OpenFF Interchange creation,
OpenMM relaxation, solvation, and downstream simulation setup. The hard part is
not only drawing a new bond. A reliable conjugation workflow must also account
for atom identity, leaving atoms, product residue definitions, charge provenance,
parameter coverage, and geometry.

The current validated path is an opt-in, config-driven NHS-Lys polymer vertical
slice. It starts from an unmodified protein, attaches an NHS-activated polymer to
lysine, ingests the product state through Pablo/OpenFF Interchange, combines
ff14SB protein treatment with polymer templates, applies a local NAGL patch
charge bridge near the linkage, reconciles local charge evidence, runs restrained
OpenMM relaxation, and finishes with final solvation.

The N-glycosylation path and multi-site attachment paths share some mechanism
and validation infrastructure. Their support level differs by scope; see
{doc}`../reference/conjugation_support_matrix` for the current matrix.

## Decision

PolyzyMD treats conjugation as a validated construction workflow rather than as
an arbitrary reaction plugin interface.

The workflow records evidence in `conjugate_validation_report.json` instead of
assuming success from a completed build. The report audits:

- product bond graph checks;
- retained atom presence and leaving atom absence;
- conservative linkage valence sanity;
- charge bridge and local reconciliation evidence;
- OpenFF Interchange to OpenMM parameter coverage;
- linkage geometry and severe close contacts;
- restrained OpenMM relaxation evidence.

This gives users and reviewers a single artifact to inspect while keeping the
implementation honest about which evidence was available, passed, warned, failed,
or was skipped.

## Consequences

The design has several consequences:

- New mechanisms cannot be supported by SMARTS alone. They must define product
  residue names, atom mappings, leaving atoms, parameterization expectations,
  charge patch behavior, and validation evidence.
- Validation is necessary but not sufficient. A passing report is an internal
  consistency audit, not a guarantee of scientific correctness for arbitrary
  chemistry.
- Pablo/OpenFF Interchange remains the product-state ingestion and
  parameterization bridge. Failures should surface as missing residue, atom,
  charge, parameterization, geometry, or relaxation evidence rather than silent
  fallbacks.
- The CUDA 12.4 simulation-only handoff remains implementation-dependent. The
  build stack and simulation runtime may need separate pixi environments until a
  stable export/reload path is documented.
- Platform fail-fast behavior is still a follow-up. Users should explicitly test
  OpenMM platform availability before production submission.

## Alternatives considered

### Allow arbitrary reaction SMARTS as a user extension point

This would make configuration look flexible, but it would hide unsupported work.
SMARTS can describe atom roles, but it does not provide coordinate placement,
product residue definitions, force-field templates, charge patches, or OpenMM
relaxation evidence. PolyzyMD therefore treats custom SMARTS as
mechanism-development input, not as a stable public plugin API.

### Rely on successful Interchange creation alone

Interchange creation is necessary, but it does not explicitly answer whether the
expected product bond is present, leaving atoms were removed, charge provenance is
acceptable, or relaxation evidence exists. The validation report keeps these
questions separate and inspectable.

### Use whole-conjugate fallback charges

Whole-conjugate AM1-BCC, Gasteiger fallback, unvalidated formal-charge
parameterization evidence, or unmarked cached charges would make builds appear
more permissive than the evidence supports. The current bridge uses mapped source
charges plus a local patch around the linkage using the pre-production OpenFF
NAGL model `openff-gnn-am1bcc-0.1.0-rc.3.pt`, then fails when required atom
identities are missing. It is not a GLYCAM, CHARMM, or AshGC workflow.

## Extension boundary

Conjugation internals are still evolving. Contributors can add mechanisms, but
should not assume a stable public mechanism-plugin API. Work with the current
internal contracts, tests, and validation checklist in
{doc}`../contributor_guide/extending_conjugation_mechanisms`, and document any
new support level in {doc}`../reference/conjugation_support_matrix`.
