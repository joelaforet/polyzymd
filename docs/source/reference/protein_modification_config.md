# Protein Modification Config Reference

This reference describes the executable protein modification configuration used
by the current `polyzymd build` workflow. Protein modifications are configured
under `conjugation:` because the supported implementation began as covalent
polymer-protein conjugation.

Current executable support covers exploratory NHS-lysine polymer attachments and
the wired moiety-plus-mechanism attachment model. Generalized recipe fields for
all protein modifications are future design material and are not part of the
current schema.

## Top-Level Block

```yaml
conjugation:
  enabled: true
  attachments: []
```

| Field | Type | Meaning |
|-------|------|---------|
| `enabled` | boolean | Turns conjugation handling on or off. This remains opt-in. |
| `attachments` | list | Ordered list of requested covalent attachments. |
| `ccd_pablo` | mapping | Pablo/CCD policy for structure ingestion and custom residues. |
| `chain_policy` | mapping | Chain assignment policy for protein and attached moieties. |
| `charge` | mapping | Charge patching and total-charge policy. |
| `diagnostics` | mapping | Controls diagnostic sidecar output. |

When `enabled: true`, at least one attachment must be present and enabled.

## Attachment

An attachment is one requested covalent residue modification. The current schema
requires an explicit `moiety` and `mechanism` for each attachment.

```yaml
attachments:
  - name: lys23-sbma-egpma-nhs
    enabled: true
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

| Field | Type | Required | Meaning |
|-------|------|----------|---------|
| `name` | string | yes | User-facing identifier for this attachment. |
| `enabled` | boolean | no | Whether this attachment should be applied. Defaults to `true`. |
| `site` | mapping | yes | Protein residue to modify. |
| `moiety` | mapping | yes | Object to attach to the protein. |
| `mechanism` | mapping | yes | Chemistry mechanism used to attach the moiety. |
| `placement` | mapping | no | Placement and clash-control settings. |

## Site

The site identifies the protein residue being modified.

```yaml
site:
  chain_id: A
  residue_name: LYS
  residue_number: 23
  atom_name: NZ
```

| Field | Type | Required | Meaning |
|-------|------|----------|---------|
| `chain_id` | string | yes | Protein chain ID. |
| `residue_number` | integer | yes | PDB residue number. |
| `residue_name` | string | recommended | Expected input residue name, such as `LYS`. |
| `atom_name` | string | mechanism-dependent | Protein atom to modify, such as `NZ`. |
| `insertion_code` | string | no | PDB insertion code. |
| `atom_serial` | integer | no | PDB atom serial for exact atom selection. |
| `atom_index` | integer | no | Zero-based atom index for exact atom selection. |

For `nhs_lys_amide`, the mechanism resolves the lysine `NZ` atom when the site
identifies a lysine residue. Include `residue_name` to make validation errors
clearer.

## Moiety

A moiety is the group added to the protein. Current executable examples use
either a polymer `recipe` or a single-residue SMILES moiety, depending on the
mechanism.

### Polymer Recipe Moiety

```yaml
moiety:
  name: SBMA-EGPMA-NHS
  recipe:
    name: SBMA-EGPMA-NHS
    length: 9
    seed: 7
    forced_reactive_monomer_label: C
    monomers:
      - label: C
        name: NHS
        residue_name: NHS
        smiles: CC(=C)C(=O)ON1C(=O)CCC1=O
        probability: 0.01
```

### SMILES Moiety

```yaml
moiety:
  name: glcnac
  smiles: CO
  residue_name: NAG
```

| Field | Type | Meaning |
|-------|------|---------|
| `name` | string | Moiety identifier. |
| `role` | string | User-facing role label, such as `polymer`, `glycan`, or `moiety`. |
| `residue_name` | string | Residue name to use for a generated single-residue moiety. |
| `input_path` | path | PDB/SDF file to load for explicit-linkage workflows. |
| `smiles` | string | SMILES for generated single-residue moieties. |
| `link_site` | mapping | Explicit moiety atom selected for an explicit linkage. |
| `recipe` | mapping | Polymer recipe for generated polymer moieties. |
| `polymer_recipe` | mapping | Equivalent explicit name for `recipe`. |

SMILES moieties used with the wired N-glycosylation mechanism must provide a
three-character PDB-safe `residue_name`.

## Mechanism

A mechanism describes how the site and moiety react.

```yaml
mechanism:
  name: nhs_lys_amide
```

| Field | Type | Meaning |
|-------|------|---------|
| `name` | string | Mechanism identifier, such as `nhs_lys_amide`, `n_glycosylation`, or `explicit_linkage`. |
| `product_residues.site` | string | Protein residue name after modification. Required for `explicit_linkage`. |
| `product_residues.moiety` | string | Moiety residue name after linkage. Required for `explicit_linkage`. |
| `bond.site_atom` | string | Protein atom participating in the new bond. Required for `explicit_linkage`. |
| `bond.moiety_atom` | string | Moiety atom participating in the new bond. Required for `explicit_linkage`. |
| `bond.order` | number | Bond order. |
| `bond.target_bond_length_angstrom` | number | Target distance after placement. |
| `leaving_atoms.site` | list | Protein atoms removed during bond formation. |
| `leaving_atoms.moiety` | list | Moiety atoms removed during bond formation. |

The `nhs_lys_amide` mechanism uses built-in NHS-lysine defaults. Explicit custom
PDB linkages must specify the atom-level fields required by the schema.

## Explicit PDB Linkage Example

Use `explicit_linkage` only when built-in mechanism defaults are insufficient.

```yaml
conjugation:
  enabled: true
  attachments:
    - name: explicit-custom-linkage
      site:
        chain_id: A
        residue_name: LYS
        residue_number: 23
        atom_name: NZ
      moiety:
        name: custom-moiety
        role: moiety
        input_path: structures/custom_moiety.pdb
        link_site:
          chain_id: X
          residue_name: LIG
          residue_number: 1
          atom_name: C1
      mechanism:
        name: explicit_linkage
        product_residues:
          site: LYX
          moiety: LIX
        bond:
          site_atom: NZ
          moiety_atom: C1
          target_bond_length_angstrom: 1.47
        leaving_atoms:
          site: [HZ1]
          moiety: [H1]
```

## Placement

```yaml
placement:
  strategy: preserve_existing
  target_bond_length_angstrom: 1.33
  clash_cutoff_angstrom: 1.5
```

| Field | Meaning |
|-------|---------|
| `strategy` | Placement strategy placeholder. |
| `target_bond_length_angstrom` | Optional target bond length for placement workflows. |
| `clash_cutoff_angstrom` | Steric clash cutoff for placement checks. |

## Support Levels

| Level | Meaning |
|-------|---------|
| Planned | Config vocabulary is documented, but execution is not implemented. |
| Validated | Config and mechanism planning work without changing topology. |
| Executable | PolyzyMD can build, relax, solvate, and write artifacts. |
| Production-ready | Templates, charges, tests, and scientific guidance are validated. |

Current support:

| Feature | Status |
|---------|--------|
| Enabled NHS-lysine polymer attachments through `moiety.recipe` and `mechanism.name: nhs_lys_amide` | Executable for exploratory workflows. |
| Multiple enabled NHS-lysine attachments | Executable for exploratory workflows. |
| SMILES moiety config with `mechanism.name: n_glycosylation` | Wired for mechanism tests and exploratory workflow development. |
| Mixed mechanisms in one config | Planned. |
| O-glycosylation | Planned. |
| General recipe field for biology-first modifications | Planned and non-executable. |

## Build Commands

Enabled conjugation configs using the current `moiety` plus `mechanism` syntax
are handled by the public conjugation workflow when you run the regular build
command:

```bash
polyzymd build -c config.yaml
```

For a conjugated GROMACS handoff, request GROMACS export from the same command:

```bash
polyzymd build -c config.yaml --format gromacs
```

The GROMACS export path completes conjugate construction, solvation, final
OpenFF Interchange creation, and then writes the `.gro`, `.top`, and `.itp`
handoff files.
