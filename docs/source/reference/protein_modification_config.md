# Protein Modification Config Reference

This reference lists the intended configuration fields for residue-level protein
modifications. It documents the target interface for the generalized framework.
The current executable implementation supports one NHS-lysine polymer attachment;
broader multi-mechanism execution is planned.

## Top-Level Block

```yaml
conjugation:
  enabled: true
  mode: construct
  attachments: []
```

| Field | Type | Meaning |
|-------|------|---------|
| `enabled` | boolean | Turns protein modification handling on or off. |
| `mode` | string | Workflow mode. `construct` builds modifications from config. |
| `attachments` | list | Ordered list of requested protein modifications. |
| `ccd_pablo` | mapping | Pablo/CCD policy for structure ingestion and custom residues. |
| `chain_policy` | mapping | Chain assignment policy for protein and attached moieties. |
| `charge` | mapping | Charge patching and total-charge policy. |
| `diagnostics` | mapping | Controls diagnostic sidecar output. |

## Attachment

An attachment is one requested residue modification.

```yaml
attachments:
  - name: lys23-polymer
    enabled: true
    site:
      chain_id: A
      residue_number: 23
    protein_modification_recipe:
      name: sbma-egpma-nhs-polymer
```

| Field | Type | Required | Meaning |
|-------|------|----------|---------|
| `name` | string | yes | User-facing identifier for this modification. |
| `enabled` | boolean | no | Whether this attachment should be applied. Defaults to `true`. |
| `site` | mapping | yes | Protein residue to modify. |
| `protein_modification_recipe` | mapping | standard workflow | Biology-level modification recipe. |
| `moiety` | mapping | advanced workflow | Explicit object to add when bypassing recipe defaults. |
| `mechanism` | mapping | advanced workflow | Explicit chemistry mechanism overrides. |
| `placement` | mapping | no | Placement and clash-control settings. |

`protein_modification_recipe` is the preferred standard interface. `moiety` and
`mechanism` remain available for advanced explicit workflows.

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
| `residue_name` | string | standard overrides | Expected input residue name, such as `LYS`. |
| `atom_name` | string | standard overrides | Protein atom to modify, such as `NZ`. |
| `insertion_code` | string | no | PDB insertion code. |

Standard recipes should infer `residue_name` and `atom_name` when the biological
request is unambiguous. Advanced users can specify them to force exact matching.

## Protein Modification Recipe

The recipe is the standard biology-first interface.

```yaml
protein_modification_recipe:
  name: phosphoserine
```

| Field | Type | Meaning |
|-------|------|---------|
| `name` | string | Recipe identifier. |
| `kind` | string | Optional recipe kind, such as `polymer`, `glycan`, `small_molecule`, or `residue_replacement`. |
| `mechanism` | string | Optional mechanism override. |
| `product_residues` | mapping | Optional output residue names. |
| `input_path` | path | Optional PDB/SDF source for the moiety. |
| `smiles` | string | Optional SMILES for generated moieties. |
| `placement` | mapping | Optional placement overrides. |
| `parameterization` | mapping | Optional force-field or template policy. |

Recipe-specific fields are allowed. Polymer recipes may include `length`,
`monomers`, probabilities, sequence controls, and Polymerist generation options.

## Moiety

A moiety is the group added to the protein.

```yaml
moiety:
  name: custom-label
  role: small_molecule
  input_path: structures/custom_label.sdf
  link_site:
    chain_id: X
    residue_name: LIG
    residue_number: 1
    atom_name: C1
```

| Field | Type | Meaning |
|-------|------|---------|
| `name` | string | Moiety identifier. |
| `role` | string | Biological role, such as `glycan`, `polymer`, `ptm`, or `moiety`. |
| `residue_name` | string | Residue name to use for the moiety when applicable. |
| `input_path` | path | PDB/SDF file to load. |
| `smiles` | string | SMILES for generated small molecules. |
| `link_site` | mapping | Explicit moiety atom selected for linkage. |
| `recipe` | mapping | Compatibility alias for polymer-style recipes. |

Prefer `protein_modification_recipe` for standard workflows. Use `moiety` when
you need explicit custom inputs.

## Mechanism

A mechanism describes how the site and moiety react.

```yaml
mechanism:
  name: nhs_lys_amide
  product_residues:
    site: LYX
    moiety: NHX
  bond:
    site_atom: NZ
    moiety_atom: C003
    order: 1
    target_bond_length_angstrom: 1.33
  leaving_atoms:
    site: [H11, H13]
    moiety: [O001, N000]
```

| Field | Type | Meaning |
|-------|------|---------|
| `name` | string | Mechanism identifier. |
| `product_residues.site` | string | Protein residue name after modification. |
| `product_residues.moiety` | string | Moiety residue name after linkage. |
| `bond.site_atom` | string | Protein atom participating in the new bond. |
| `bond.moiety_atom` | string | Moiety atom participating in the new bond. |
| `bond.order` | number | Bond order. |
| `bond.target_bond_length_angstrom` | number | Target distance after placement. |
| `leaving_atoms.site` | list | Protein atoms removed during bond formation. |
| `leaving_atoms.moiety` | list | Moiety atoms removed during bond formation. |

Standard recipes should fill these fields automatically.

## Placement

```yaml
placement:
  reactive_sphere_radius_angstrom: 2.0
  tolerance_angstrom: 2.0
  movebadrandom: true
  nloop: 1000
```

| Field | Meaning |
|-------|---------|
| `reactive_sphere_radius_angstrom` | Radius around the protein site where the moiety reactive atom is placed. |
| `tolerance_angstrom` | Packmol distance tolerance. |
| `target_bond_length_angstrom` | Final snapped bond length. |
| `movebadrandom` | Allows Packmol to move badly packed structures randomly. |
| `nloop` | Maximum Packmol optimization loops. |

Large moieties should use placement retries and clash scoring before OpenMM
relaxation.

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
| One NHS-lysine polymer attachment | Executable for exploratory workflows. |
| `protein_modification_recipe` standard field | Planned. |
| Multiple attachments | Planned. |
| Mixed mechanisms in one config | Planned. |
| O-glycosylation | Planned. |
| N-glycosylation | Planned mechanism placeholder. |
| Phosphorylation | Planned. |
| PEGylation | Planned. |

## Standard Example

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

    - name: asn67-n-glycan
      site:
        chain_id: A
        residue_number: 67
      protein_modification_recipe:
        name: high-mannose-n-glycosylation
```

## Advanced Explicit Example

```yaml
conjugation:
  enabled: true
  mode: construct
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

Use explicit configuration only when recipe defaults are insufficient.
