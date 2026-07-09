# Protein Modification Architecture

PolyzyMD's protein modification framework is intended to let users describe
biology in a configuration file and let the software handle the molecular
modeling details. A user should be able to say, "add this sugar to this residue"
or "attach this polymer to this lysine" without first becoming an expert in atom
names, graph edits, force-field templates, or coordinate generation.

The current executable implementation supports one NHS-lysine polymer
conjugation. The architecture described here is the reference design for
generalizing that implementation to arbitrary residue modifications, including
glycosylation, phosphorylation, PEGylation, polymer attachment, lipidation, and
custom covalent labels.

## Design Goal

The future standard workflow should be biology-first. The example below is
design material and is not accepted by the current executable schema:

```yaml
conjugation:
  enabled: true
  attachments:
    - name: ser4-o-glycan
      site:
        chain_id: A
        residue_number: 4
      protein_modification_recipe:
        name: core1-o-glycosylation
```

In that future standard workflow, PolyzyMD should resolve the residue identity, choose
the correct atom, load or generate the modification, apply the mechanism, place
the added group, relax the structure, and write diagnostics explaining what it
did.

The current executable workflow uses the same top-level `conjugation.enabled`
and `conjugation.attachments` block, but each attachment must provide explicit
`moiety` and `mechanism` fields.

Advanced users should still be able to override the resolved atoms, provide
custom PDB or SDF inputs, specify product residue names, define leaving atoms,
and supply custom mechanism or parameterization settings.

## Vocabulary

These words are used throughout the protein modification workflow.

| Word | Meaning | Example |
|------|---------|---------|
| Protein modification | Any covalent change to a protein residue. | Glycosylation, phosphorylation, PEGylation, polymer attachment. |
| Attachment | One requested modification in the config file. | `lys23-polymer` or `ser4-o-glycan`. |
| Site | The protein residue being modified. | Chain `A`, residue `LYS 23`. |
| Moiety | The chemical group being added to the protein. | A sugar, phosphate, PEG chain, polymer, fluorophore, or lipid. |
| Mechanism | The chemistry rule for how the moiety attaches. | NHS ester plus lysine forms an amide. |
| Protein modification recipe | A reusable recipe that knows what to add and which mechanism to use. | `core1-o-glycosylation` or `sbma-egpma-nhs-polymer`. |
| Product residue | The residue name after modification. | Phosphoserine may become `SEP`; linked lysine may become `LYX`. |
| Leaving atoms | Atoms removed during bond formation. | A lysine hydrogen or an NHS leaving group. |
| Placement | Coordinate generation for the added moiety near the target site. | Packmol places a polymer near lysine `NZ`. |
| Relaxation | Energy minimization or short MD after placement. | OpenMM minimization after a polymer is linked. |

## Why Plan `protein_modification_recipe`

A future recipe field would package the biological intent into a reusable
object. Instead of forcing every user to define atom-level chemistry, a recipe
could provide defaults for the moiety, mechanism, product residue names,
placement policy, and parameterization expectations.

This planned, non-executable syntax would keep common workflows short:

```yaml
protein_modification_recipe:
  name: phosphoserine
```

It would also keep advanced workflows possible:

```yaml
protein_modification_recipe:
  name: custom-nhs-polymer
  kind: polymer
  length: 31
  mechanism: nhs_lys_amide
  product_residues:
    site: LYX
    moiety: NHX
  placement:
    reactive_sphere_radius_angstrom: 2.0
```

## Conceptual Data Flow

The generalized workflow should follow this sequence:

```text
config.yaml
  -> attachment list
  -> protein modification recipe registry
  -> moiety builder
  -> reaction template
  -> atom-level attachment plan
  -> placement
  -> graph or PDB assembly
  -> parameterization
  -> relaxation
  -> final modified topology/PDB
  -> solvation and optional free-polymer packing
```

Each step should write enough diagnostics that a non-expert can answer four
questions:

1. Which protein residue was modified?
2. Which moiety was added?
3. Which bond was created?
4. What assumptions were made about charges, atom names, and templates?

## Main Components

### Attachment Planner

The planner reads `conjugation.attachments` and turns each enabled attachment
into an ordered work item. It should validate that each requested site exists,
that the mechanism is compatible with the residue, and that two attachments do
not accidentally target the same residue unless explicitly allowed.

### Recipe Registry

The recipe registry maps user-facing names to modification defaults. A recipe is
not a chemistry executor by itself. It answers questions like:

| Question | Example answer |
|----------|----------------|
| What kind of moiety is this? | `polymer`, `glycan`, `small_molecule`, `residue_replacement`. |
| How should the moiety be built? | Generate from monomers, load PDB, load SDF, or use a template. |
| Which mechanism attaches it? | `nhs_lys_amide`, `o_glycosidic_ser_thr`, `phosphorylation`. |
| What defaults are safe? | Product residues, bond length, placement settings, charge policy. |

### Moiety Builder

The moiety builder creates or loads the object that will be attached. Current
polymer work uses Polymerist to generate polymer fragments. Future builders may
load glycans from curated PDB/SDF files, create PEG chains, or replace one amino
acid residue with a modified residue template.

### Reaction Templates

Reaction templates define how a class of chemistry works. A template should
declare:

| Field | Purpose |
|-------|---------|
| Allowed sites | Which protein residue and atom can react. |
| Moiety reactive group | Which atom on the added group forms the new bond. |
| Leaving atoms | Which atoms are removed from protein and moiety. |
| Product residue names | How the final linked residues should be named. |
| Bond specification | Bond order and target bond length. |
| Parameterization policy | Whether curated templates are required. |

### Placement Engine

The placement engine generates starting coordinates for the added moiety. The
current workflow uses Packmol followed by OpenMM minimization. Larger moieties
need robust placement scoring: multiple Packmol seeds, clash checks after bond
snapping, and clear failure messages when no acceptable placement is found.

### Parameterization Layer

The parameterization layer decides whether a result is suitable for exploratory validation
or production simulations. Some modifications can use generic OpenFF-style
fallbacks for exploratory runs. Others, especially glycans and post-translational
modifications with known residue templates, should require curated templates and
validated charges.

## Multiple Attachments

The schema already represents multiple attachments as a list. The generalized
executor should apply them sequentially:

```text
input protein
  -> attachment 1 product
  -> attachment 2 product
  -> attachment 3 product
  -> relaxed multi-modified protein
```

Sequential execution matters because each modification changes the structure the
next modification sees. For example, if a protein is O-glycosylated at residue 4,
then polymer-modified at residue 23, and then N-glycosylated at residue 67, the
second and third placement steps should avoid clashes with earlier
modifications.

## Support Levels

Documentation and diagnostics should use explicit support levels.

| Level | Meaning |
|-------|---------|
| Planned | The config can express the workflow, but execution is not implemented. |
| Validated | PolyzyMD can parse and validate the attachment plan without changing topology. |
| Executable | PolyzyMD can build, relax, solvate, and write artifacts for the workflow. |
| Production-ready | The workflow has validated templates, charges, tests, and scientific guidance. |

Current status:

| Workflow | Status |
|----------|--------|
| One NHS-lysine polymer attachment | Executable for exploratory workflows. |
| Multiple NHS-lysine attachments | Planned. |
| O-glycosylation | Planned. |
| N-glycosylation | Planned mechanism placeholder exists. |
| Phosphorylation | Planned. |
| PEGylation | Planned. |
| Mixed mechanisms in one config | Planned. |

## Design Principles

1. Prefer biology-level configuration for standard workflows.
2. Preserve atom names and residue identities whenever possible.
3. Record every resolved atom-level decision in diagnostics.
4. Separate moiety generation from reaction mechanism logic.
5. Separate placement from relaxation.
6. Treat generic charge fallbacks as exploratory unless validated templates are available.
7. Keep advanced overrides available without making them mandatory for common cases.

## Implementation Implications

Generalization should not add a new one-off path for every biology workflow.
Instead, new capabilities should enter through the recipe registry, mechanism
registry, and moiety builders. The NHS-lysine polymer path should become the
first executable specialization of the same framework used by glycosylation,
phosphorylation, PEGylation, and other residue modifications.
