# Extend Protein Modifications

Use this guide when you want to add a new kind of protein residue modification
to PolyzyMD. The goal is to expose a simple recipe name to standard users while
keeping atom-level chemistry available for advanced users.

This is an implementation guide for the generalized framework. Some interfaces
described here are planned rather than fully executable today.

## Decide The User Story First

Start from the biology sentence a user should be able to write:

```text
Add <modification> to <residue> of my protein.
```

Examples:

| User story | Config recipe name |
|------------|--------------------|
| Add phosphate to serine 12. | `phosphoserine` |
| Add an O-linked core 1 glycan to serine 4. | `core1-o-glycosylation` |
| Attach a PEG linker to cysteine 20. | `cys-pegylation` |
| Attach an NHS-reactive polymer to lysine 23. | `sbma-egpma-nhs-polymer` |

The standard config should not require users to know the reactive atom unless
the residue is ambiguous or the workflow is custom.

## Choose A Support Level

Before adding code, decide what level you are implementing.

| Level | Required behavior |
|-------|-------------------|
| Planned | Document the intended config and scientific assumptions. |
| Validated | Parse config and validate sites, recipes, and mechanisms. |
| Executable | Build coordinates, create bonds, relax, and write artifacts. |
| Production-ready | Provide validated templates, charges, tests, and usage guidance. |

Do not present a workflow as production-ready if it relies on exploratory charge
fallbacks or unvalidated residue templates.

## Add A Protein Modification Recipe

A `protein_modification_recipe` should package defaults for standard users. It
should answer:

| Question | Example |
|----------|---------|
| What is added? | phosphate, glycan, PEG, polymer, label |
| How is it built? | built-in template, PDB input, SDF input, polymer generation |
| How does it attach? | mechanism name |
| Which residues are allowed? | serine/threonine, lysine, asparagine |
| What product names are used? | `SEP`, `LYX`, `NGX` |
| What parameterization is required? | curated template, OpenFF fallback, custom charges |

The target config shape is:

```yaml
attachments:
  - name: ser12-phosphorylation
    site:
      chain_id: A
      residue_number: 12
    protein_modification_recipe:
      name: phosphoserine
```

Only expose low-level fields when users need to override recipe defaults.

## Add Or Reuse A Reaction Template

A reaction template describes how the modification attaches. Add a new template
when a new reaction class is needed. Reuse an existing template when only the
moiety changes.

A reaction template should define:

| Template field | Purpose |
|----------------|---------|
| `identifier` | Stable name used in config and diagnostics. |
| `allowed_sites` | Residue and atom rules for the protein site. |
| `moiety_reactive_group` | Reactive atom or group on the moiety. |
| `graph_edits` | Atoms removed and bonds added. |
| `charge_patch_hint` | Expected charge/template strategy. |
| `rationale` | Plain-language scientific explanation. |

For built-ins, add a `ReactionTemplate` implementation under `reactions/` and
include tests that it exposes stable defaults and rejects incompatible inputs.

## Add A Moiety Builder If Needed

Use the simplest builder that preserves chemical intent.

| Moiety source | Recommended builder |
|---------------|---------------------|
| Curated modified residue | Residue template or PDB/SDF loader. |
| Glycan | Curated glycan PDB/SDF loader first, recipe generator later. |
| Polymer | Polymerist-backed polymer generator. |
| PEG or linker | Small-molecule SDF loader or deterministic generator. |
| Fluorophore or custom label | User-supplied SDF/PDB with explicit link site. |

The moiety builder should produce a structure with atom names stable enough for
mechanism resolution and diagnostics.

## Define Product Residues And Atom Names

Product residue names are part of the modeling contract. They affect PDB output,
Pablo/CCD ingestion, force-field matching, visualization, and analysis.

Use standard residue names when they exist, such as `SEP` for phosphoserine.
Use custom names only when no standard template exists, and document what they
mean.

When possible, preserve input protein atom names in output PDB files so tools
like PyMOL can still identify the protein backbone and render cartoon
representations.

## Add Placement Policy

Large moieties need placement safeguards. A robust executable workflow should:

1. Try multiple Packmol seeds or placements.
2. Snap the reactive atoms to the target bond length.
3. Score the snapped structure for protein-moiety clashes.
4. Retry or fail before OpenMM if severe clashes remain.
5. Run OpenMM minimization after placement.

Do not rely on Polymerist or Open Babel minimization for large conjugated
polymers when the intended relaxation path is OpenMM.

## Add Parameterization Policy

Every executable mechanism needs a parameterization decision.

| Policy | Appropriate use |
|--------|-----------------|
| Curated residue template | Standard PTMs, glycans, production workflows. |
| Custom force-field fragment | New residue names or custom linkers. |
| Direct OpenFF fallback | Exploratory smoke tests with clear warnings. |
| User-supplied charges | Advanced workflows with known charge models. |

Diagnostics should explicitly state whether the build is exploratory or
production-ready.

## Add Tests

At minimum, add tests for:

| Test type | What it should prove |
|-----------|----------------------|
| Config parsing | The recipe and attachment block validate. |
| Mechanism validation | Allowed sites pass and incompatible sites fail. |
| Moiety resolution | The reactive atom and leaving atoms are resolved. |
| Placement | Severe clashes are detected or retried. |
| PDB output | Protein atom names and chain IDs are preserved. |
| Diagnostics | Assumptions and support level are reported. |

For multiple-attachment workflows, include a sequential test where attachment 2
uses the product of attachment 1 as its input structure.

## Document The Standard Workflow

When you add a new recipe, update user-facing docs before calling it supported.
The docs should include:

- a plain-language description of the modification
- compatible residue types
- minimal config
- expected output residue names
- limitations and support level
- any template or charge requirements

Keep the standard workflow short. Put atom-level overrides in reference or
advanced how-to sections, not in the first example.

## Checklist

Before marking a modification executable, confirm:

- standard config works without atom-level overrides
- advanced config can override link atoms and leaving atoms
- diagnostics explain resolved atoms and created bonds
- output PDB preserves protein atom names
- placement has clash checks or retries
- OpenMM relaxation succeeds on a small fixture
- parameterization limitations are explicit
- docs state the correct support level
