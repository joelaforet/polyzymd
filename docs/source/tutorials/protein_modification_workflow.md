# Build a Protein With Multiple Modifications

This tutorial shows the intended user experience for building a protein with
multiple residue modifications in one `config.yaml`. It is written as a
reference workflow for implementation: the config shape is the target interface,
while full multi-mechanism execution is still planned.

By the end, you will have a config that describes three biological requests:

1. Add an O-linked glycan to residue 4.
2. Attach a polymer to residue 23.
3. Add an N-linked glycan to residue 67.

## Before You Start

You need:

- a prepared protein PDB
- residue numbers for the modifications you want
- names of the modifications you want to apply
- optional PDB or SDF files for custom moieties

You do not need to know atom names for standard built-in recipes. The recipe and
mechanism registries should resolve those for you.

```{important}
Current executable support is narrower than this tutorial. Today, PolyzyMD can
execute one NHS-lysine polymer attachment through the prototype workflow. This
page defines the intended multi-modification workflow that the implementation
should grow into.
```

## Step 1: Start With The Protein

Begin with the normal PolyzyMD protein block:

```yaml
name: multi-modified-protein

enzyme:
  name: example-protein
  pdb_path: structures/protein.pdb
```

The PDB should already have the residue numbering you want to use. For example,
if you say `residue_number: 4`, PolyzyMD will look for residue 4 in the selected
chain.

## Step 2: Enable Protein Modifications

Protein modifications live under `conjugation:` because the original executable
path began as polymer-protein conjugation. The generalized framework treats this
section as the place for any covalent residue modification.

```yaml
conjugation:
  enabled: true
  mode: construct
  attachments: []
```

`mode: construct` means PolyzyMD should build the modified structure from the
unmodified input protein and the requested attachments.

## Step 3: Add O-Glycosylation At Residue 4

For a standard recipe, the user should only need the site and recipe name:

```yaml
conjugation:
  enabled: true
  mode: construct
  attachments:
    - name: ser4-o-glycan
      site:
        chain_id: A
        residue_number: 4
      protein_modification_recipe:
        name: core1-o-glycosylation
```

The recipe should provide the default mechanism, allowed residue types, moiety
structure, and product residue naming. If residue 4 is not a compatible residue,
validation should fail with a clear message before any modeling work begins.

## Step 4: Add A Polymer At Residue 23

Polymer attachment can use the same standard recipe pattern. The recipe can
include polymer-specific options such as length and sequence controls:

```yaml
    - name: lys23-polymer
      site:
        chain_id: A
        residue_number: 23
      protein_modification_recipe:
        name: sbma-egpma-nhs-polymer
        length: 31
        forced_reactive_monomer_label: C
        forced_reactive_monomer_index: 15
        monomers:
          - label: A
            name: SBMA
            residue_name: SBM
            probability: 0.945
            smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
          - label: B
            name: EGPMA
            residue_name: EGP
            probability: 0.045
            smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
          - label: C
            name: NHS
            residue_name: NHS
            probability: 0.01
            smiles: "CC(=C)C(=O)ON1C(=O)CCC1=O"
```

The NHS monomer is the reactive monomer. The mechanism should know that an NHS
ester can react with lysine and should resolve the lysine `NZ` atom and the NHS
acyl carbon.

## Step 5: Add N-Glycosylation At Residue 67

Add the third attachment to the same list:

```yaml
    - name: asn67-n-glycan
      site:
        chain_id: A
        residue_number: 67
      protein_modification_recipe:
        name: high-mannose-n-glycosylation
```

For N-glycosylation, the mechanism should validate that residue 67 is a
compatible asparagine site. A future implementation may also validate the
sequence context around the asparagine when that is scientifically required.

## Step 6: Review The Complete Block

The complete modification block is:

```yaml
conjugation:
  enabled: true
  mode: construct
  attachments:
    - name: ser4-o-glycan
      site:
        chain_id: A
        residue_number: 4
      protein_modification_recipe:
        name: core1-o-glycosylation

    - name: lys23-polymer
      site:
        chain_id: A
        residue_number: 23
      protein_modification_recipe:
        name: sbma-egpma-nhs-polymer
        length: 31
        forced_reactive_monomer_label: C
        forced_reactive_monomer_index: 15
        monomers:
          - label: A
            name: SBMA
            residue_name: SBM
            probability: 0.945
            smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
          - label: B
            name: EGPMA
            residue_name: EGP
            probability: 0.045
            smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
          - label: C
            name: NHS
            residue_name: NHS
            probability: 0.01
            smiles: "CC(=C)C(=O)ON1C(=O)CCC1=O"

    - name: asn67-n-glycan
      site:
        chain_id: A
        residue_number: 67
      protein_modification_recipe:
        name: high-mannose-n-glycosylation
```

Attachments are applied in order. This means the polymer placement should see
the O-glycan, and the N-glycan placement should see both earlier modifications.

## Step 7: Add Solvent And Optional Free Polymers

Protein modifications define covalent changes to the protein. Free polymers are
still separate non-covalent solutes and should stay under `polymers:`.

```yaml
polymers:
  enabled: true
  generation_mode: dynamic
  length: 3
  count: 1

solvent:
  primary:
    model: tip3p
  ions:
    neutralize: true
  box:
    shape: cube
    padding: 0.5
```

This distinction is important. A polymer in `conjugation.attachments` is bonded
to the protein. A polymer in `polymers:` is packed near the protein but remains
unattached.

## Step 8: Validate Before Building

The generalized workflow should provide a validation command before expensive
coordinate generation:

```bash
polyzymd validate -c config.yaml
```

Validation should report:

| Check | Expected result |
|-------|-----------------|
| Site exists | Each requested residue is present in the PDB. |
| Mechanism compatibility | Each residue can use the requested mechanism. |
| Recipe resolution | Each recipe name maps to a known modification recipe. |
| Attachment conflicts | No two attachments unexpectedly modify the same residue. |
| Template readiness | Required force-field or residue templates are available. |

## Step 9: Build And Inspect Diagnostics

The eventual build command should write a modified protein PDB, a relaxed PDB,
and diagnostic JSON files:

```bash
polyzymd build -c config.yaml
```

Diagnostics should answer these questions in plain language:

1. Which residues were modified?
2. Which atoms were bonded?
3. Which atoms were removed?
4. Which residue names changed?
5. Which templates or charge policies were used?
6. Which assumptions are exploratory rather than production-ready?

## Expected Output

A successful multi-modification build should produce artifacts like:

```text
artifacts/multi-modified-protein/
|- conjugation_diagnostics.json
|- protein_modification_workflow.json
|- attachment-001-ser4-o-glycan/
|- attachment-002-lys23-polymer/
|- attachment-003-asn67-n-glycan/
|- relaxed_modified_protein.pdb
`- solvated_modified_system.pdb
```

Each attachment directory should contain the local placement inputs, linked PDB,
parameterization report, and relaxation report for that modification.

## What Works Today

The current executable workflow supports a narrower path:

| Requested workflow | Current status |
|--------------------|----------------|
| One NHS-lysine polymer attachment | Executable for exploratory builds. |
| Multiple attachments in one config | Planned. |
| O-glycosylation recipe | Planned. |
| N-glycosylation recipe | Planned. |
| Mixed glycan plus polymer mechanisms | Planned. |

Use this tutorial as the target user experience when extending the framework.
For the currently executable NHS-lysine polymer path, see the conjugation
notebook and generated artifacts used during development.
