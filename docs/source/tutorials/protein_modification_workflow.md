# Build a Protein With an NHS-Lysine Polymer Attachment

This tutorial shows the current validated vertical slice for protein
modification: an NHS-reactive polymer attached to a lysine through
`conjugation.attachments`. The examples use the schema accepted by
`polyzymd build` today.

By the end, you will have a `config.yaml` fragment that:

- enables the conjugation workflow
- selects a lysine residue on chain `A`
- defines the polymer as a `moiety.recipe`
- requests the `nhs_lys_amide` mechanism
- can be used with the regular build and GROMACS export commands

```{important}
General biology-first recipe syntax for arbitrary protein modifications is
planned but not executable in the current schema. This tutorial intentionally
uses the supported `moiety` plus `mechanism` syntax.
```

## Before You Start

You need:

- a prepared protein PDB with residue numbering you trust
- the chain ID and residue number for a lysine modification site
- a polymer recipe with an NHS-containing reactive monomer

The examples below use chain `A`, residue `23`, and an SBMA/EGPMA/NHS polymer.

## Step 1: Start With the Protein

Begin with the normal PolyzyMD protein block:

```yaml
name: lysine-polymer-conjugate

enzyme:
  name: example-protein
  pdb_path: structures/protein.pdb
```

The PDB should already use the residue numbering you want. If you set
`residue_number: 23`, PolyzyMD looks for residue 23 on the selected chain.

## Step 2: Enable Conjugation

Protein modifications live under `conjugation:` because the executable path is
the conjugation workflow.

```yaml
conjugation:
  enabled: true
  attachments: []
```

When this block contains enabled attachments,
`pixi run -e build polyzymd build -c config.yaml` routes the build through
conjugate construction before solvation and export.

## Step 3: Select the Lysine Site

Add one attachment and identify the target residue:

```yaml
conjugation:
  enabled: true
  attachments:
    - name: lys23-sbma-egpma-nhs
      site:
        chain_id: A
        residue_name: LYS
        residue_number: 23
```

The `nhs_lys_amide` mechanism resolves the lysine `NZ` atom from this site. The
`residue_name` field is recommended because it catches accidental residue-number
mismatches early.

## Step 4: Add the Polymer Moiety

Define the polymer under `moiety.recipe`. The NHS monomer is labeled `C` and is
forced into the generated sequence so the linker has one reactive site.

```yaml
      moiety:
        name: SBMA-EGPMA-NHS
        recipe:
          name: SBMA-EGPMA-NHS
          length: 9
          seed: 7
          forced_reactive_monomer_label: C
          monomers:
            - label: A
              name: SBMA
              residue_name: SBM
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
              probability: 0.945
            - label: B
              name: EGPMA
              residue_name: EGP
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
              probability: 0.045
            - label: C
              name: NHS
              residue_name: NHS
              smiles: CC(=C)C(=O)ON1C(=O)CCC1=O
              probability: 0.01
```

You can use `polymer_recipe:` instead of `recipe:` if you prefer the explicit
field name. Both keys map to the same current schema field.

## Step 5: Choose the Mechanism

Set the mechanism name to the built-in NHS-lysine amide linkage:

```yaml
      mechanism:
        name: nhs_lys_amide
```

This mechanism applies the current NHS-lysine defaults for the reactive lysine
atom, NHS reactive group, leaving atoms, product residue names, and target bond
length.

## Step 6: Review the Complete Block

The complete executable conjugation block is:

```yaml
conjugation:
  enabled: true
  attachments:
    - name: lys23-sbma-egpma-nhs
      site:
        chain_id: A
        residue_name: LYS
        residue_number: 23
      moiety:
        name: SBMA-EGPMA-NHS
        recipe:
          name: SBMA-EGPMA-NHS
          length: 9
          seed: 7
          forced_reactive_monomer_label: C
          monomers:
            - label: A
              name: SBMA
              residue_name: SBM
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
              probability: 0.945
            - label: B
              name: EGPMA
              residue_name: EGP
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
              probability: 0.045
            - label: C
              name: NHS
              residue_name: NHS
              smiles: CC(=C)C(=O)ON1C(=O)CCC1=O
              probability: 0.01
      mechanism:
        name: nhs_lys_amide
```

Multiple NHS-lysine polymer attachments can be listed in order when each entry
uses the same current attachment shape.

## Step 7: Add Solvent and Optional Free Polymers

Covalent polymer attachments belong under `conjugation.attachments`. Free
polymers are separate non-covalent solutes and stay under `polymers:`.

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
    shape: orthorhombic
    padding: 0.5
```

A polymer in `conjugation.attachments` is bonded to the protein. A polymer in
`polymers:` is packed near the protein but remains unattached.

## Step 8: Build and Inspect Diagnostics

Run the regular build command from your project root:

```bash
pixi run -e build polyzymd build -c config.yaml
```

For a conjugated GROMACS handoff, request GROMACS export from the same build
command:

```bash
pixi run -e build polyzymd build -c config.yaml --format gromacs
```

The GROMACS path first completes the conjugation workflow, then exports the final
OpenFF Interchange to `.gro`, `.top`, and `.itp` files.

Diagnostics should answer these questions in plain language:

1. Which residue was modified?
2. Which polymer moiety was added?
3. Which atoms were bonded?
4. Which atoms were removed?
5. Which residue names changed?
6. Which validation report sections passed, warned, failed, or were skipped?

## Expected Output

A successful build writes artifacts similar to:

```text
artifacts/lysine-polymer-conjugate/
|- conjugate-construction/
|  |- assembled_crosslinked.pdb
|  |- conjugate_relaxation.json
|  |- conjugate_relaxed.pdb
|  `- conjugate_validation_report.json
|- conjugated_polymer_system_workflow.json
|- system.xml
`- solvated_conjugate_free_polymers.pdb
```

Exact filenames can vary by workflow branch and export options, but the result
object records the crosslinked, relaxed, solvated, validation, and export paths
that were created.

## What Works Today

| Requested workflow | Current status |
|--------------------|----------------|
| Enabled NHS-lysine polymer attachments with `moiety.recipe` | Validated vertical slice for the stated NHS-Lys polymer scope. |
| Multiple enabled NHS-lysine polymer attachments in one config | Executable; validation metadata is multi-site aware, but scientific validation remains system-specific. |
| General biology-first recipe field | Planned and non-executable. |
| O-glycosylation recipe | Planned. |
| Mixed glycan plus polymer mechanisms | Planned. |

Use the current `moiety` plus `mechanism` attachment shape for examples paired
with `polyzymd build`.
