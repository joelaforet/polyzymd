# Protein Modification Config Reference

This reference describes the executable protein modification configuration used
by the current `polyzymd build` workflow. Protein modifications are configured
under `conjugation:` because the supported implementation began as covalent
polymer-protein conjugation.

Current validated vertical-slice support covers NHS-lysine polymer attachments.
The active architecture is generic: a moiety provider resolves each source, a
mechanism resolves a product-state attachment plan, and the workflow carries an
attachment spec through shared assembly, Pablo/OpenFF parameterization, charge
patching, conjugate relaxation, and validation. Non-NHS mechanisms and generic explicit
linkages are executable development paths, but they remain experimental unless a
mechanism-specific validation note says otherwise. Phase 11 also supports a
coordinate-only GlyGen/GlyCAM-style PDB N-glycan input path for
`mechanism.name: n_glycosylation`; its Pablo/OpenFF continuation is explicitly
experimental.

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

Internally, each enabled attachment becomes an attachment spec containing the
resolved moiety source, resolved attachment plan, product residue mapping, and
sidecar paths. Multi-attachment builds assemble all specs into one product PDB
before relaxation and validation run, so shared validation and charge-patch logic
sees the combined conjugate.

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

A moiety is the group added to the protein. The active provider layer resolves
polymer recipes, single-residue SMILES moieties, and selected file-based sources.
In Phase 11, `moiety.input_path` is used for GlyGen/GlyCAM-style residue-resolved
PDB N-glycan input with `mechanism.name: n_glycosylation`; it is also part of the
advanced explicit-linkage schema. The mechanism decides which source forms are
supported for a given attachment.

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

### GlyGen/GlyCAM PDB N-glycan Moiety

For the Phase 11 GlyGen PDB N-glycosylation path, provide the glycan PDB with
`moiety.input_path` and use `mechanism.name: n_glycosylation`. Do not also set
`smiles`, `residue_name`, `recipe`, or `polymer_recipe` on the same moiety;
the provider requires exactly one source. The loaded PDB must contain one
connected residue-resolved glycan graph and the strict reducing-end `ROH`
leaving group with atoms named `O1` and `HO1`.

```yaml
moiety:
  name: G80966KZ
  role: glycan
  input_path: structures/G80966KZ_glycam.pdb
```

The default workflow mode writes a coordinate-only artifact and a GlyGen
ingestion sidecar. It does not produce a final OpenMM system, GLYCAM/CHARMM
parameters, or a production glycan force-field assignment. See
{doc}`../how_to/glygen_glycan_conjugation` for the task workflow.

| Field | Type | Meaning |
|-------|------|---------|
| `name` | string | Moiety identifier. |
| `role` | string | User-facing role label, such as `polymer`, `glycan`, or `moiety`. |
| `residue_name` | string | Residue name to use for a generated single-residue moiety. |
| `input_path` | path | PDB file for Phase 11 GlyGen/GlyCAM-style N-glycan input with `mechanism.name: n_glycosylation`; also used by the advanced explicit-linkage schema. |
| `smiles` | string | SMILES for generated single-residue moieties. |
| `link_site` | mapping | Explicit moiety atom selected for an explicit linkage. |
| `recipe` | mapping | Polymer recipe for generated polymer moieties. |
| `polymer_recipe` | mapping | Equivalent explicit name for `recipe`. |

SMILES moieties used with the experimental N-glycosylation mechanism must
provide a three-character PDB-safe `residue_name`.

GlyGen/GlyCAM PDB N-glycan moieties used with `n_glycosylation` derive their
reactive selector from the loaded glycan: chain `C`, the reducing sugar residue,
atom `C1`, and the source atom serial. The workflow removes the glycan
`ROH:O1`/`ROH:HO1` leaving atoms and links the glycan reducing-end `C1` to the
Asn site atom.

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

The `nhs_lys_amide` mechanism uses built-in NHS-lysine defaults and is the
validated vertical slice. Experimental non-NHS mechanisms use the same resolved
plan and attachment-spec pipeline, but users should validate product chemistry,
charges, and geometry carefully. Explicit custom PDB linkages must specify the
atom-level fields required by the schema.

The public engine also accepts direct requests built from `protein_pdb_path` plus
the same `attachments` entries. Raw molecule/topology keyword arguments are not
supported by the public API.

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

## Product-State Charge Bridge

Final OpenFF Interchange creation for a covalent protein-polymer conjugate uses
an explicit full-molecule charge template. OpenFF sees the conjugated protein and
attached polymer as one covalent molecule, so PolyzyMD must pass one complete
`charge_from_molecules` template for that molecule rather than separate protein
and polymer fragments.

The current bridge is an active product-state interoperability bridge for the
validated NHS-lysine vertical slice, not a whole-conjugate AM1-BCC model and not
a residue-resolved GLYCAM or CHARMM conjugate-template workflow. Its model choices are:

- Standard protein atoms use ff14SB-style charges from the prepared source
  protein.
- Free polymers, substrate, water, and co-solvents continue to use the existing
  standard charged-template path.
- Attached polymer or moiety atoms are included in a private product patch when
  they are part of the modified residue context. The final conjugate template is
  accepted only when every atom has explicit charge provenance.
- Modified-residue and moiety atoms are charged through a private peptide-capped
  product patch: ACE--GLY--modified residue with the attached moiety--GLY--NME.
  PolyzyMD builds this bounded patch internally, charges it with the bundled
  OpenFF NAGL model `openff-gnn-am1bcc-0.1.0-rc.3.pt`, and maps charges back
  only to real product atoms. This NAGL model is a pre-production default; it is
  not GLYCAM, CHARMM, or AshGC, and the bridge does not claim those validation
  scopes.
- Any residual total-charge closure is applied only over real mapped local patch
  atoms and only within fixed internal safety bounds. Users cannot relax or tune
  the reconciliation threshold from configuration, environment variables, or API
  arguments.

The final production Interchange path refuses whole-conjugate AM1-BCC,
Gasteiger fallback, charge templates without validation evidence, and unmarked cached charges.
If any required atom cannot be mapped to one of the source classes above, the
build fails with an explicit missing-identity message rather than inventing a
fallback charge.

The bridge writes `conjugate-construction/product_state_charge_bridge.json`.
Inspect this sidecar to confirm the pre-production NAGL model, ff14SB atom
count, local patch atom count, total charge, formal charge, and any
normalization correction.

The local patch builder uses `ResolvedAttachmentPlan`, attachment specs, the
Pablo crosslink requirement, product residue mappings, generated-fragment atom and
bond-order metadata, and leaving-atom metadata. It removes leaving atoms, builds
the fixed peptide-capped ACE--GLY--modified-residue(moiety)--GLY--NME product
patch, charges that local product-state molecule with OpenFF NAGL
`openff-gnn-am1bcc-0.1.0-rc.3.pt`, and maps charges back only to real product
atoms. This is not an AshGC, GLYCAM, or CHARMM charge assignment. If the mechanism
metadata is not specific enough to build the product patch or if final topology
atom identities are missing, PolyzyMD fails clearly instead of falling back to
hardcoded NHS-lysine atom names, raw sidecar charges, or count-only ordered charge
transfer.

## Support Levels and Validation Reports

See {doc}`conjugation_support_matrix` for the current support matrix. In brief,
the validated reliability milestone covers the opt-in config-driven and direct
request NHS-lysine polymer vertical slice, including Pablo/OpenFF Interchange,
ff14SB plus polymer templates, local pre-production NAGL patch charge bridge,
local charge reconciliation, conjugate relaxation, final solvation, and
`conjugate_validation_report.json`.

The validation report audits product bond graph evidence, link atom presence,
leaving atom absence, valence sanity, charge evidence, parameter coverage,
linkage geometry, OpenMM relaxation evidence, and paths to supporting sidecars.
Validation is necessary but not sufficient: a passing report is an internal
consistency audit, not a guarantee for arbitrary chemistry.

## Build Commands

Enabled conjugation configs using the current `moiety` plus `mechanism` syntax
are handled by the public conjugation workflow when you run the regular build
command:

```bash
pixi run -e build polyzymd build -c config.yaml
```

For a conjugated GROMACS handoff, request GROMACS export from the same command:

```bash
pixi run -e build polyzymd build -c config.yaml --format gromacs
```

The GROMACS export path completes conjugate construction, solvation, final
OpenFF Interchange creation, and then writes the `.gro`, `.top`, and `.itp`
handoff files.
