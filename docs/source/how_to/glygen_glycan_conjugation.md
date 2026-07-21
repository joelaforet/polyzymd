# Build an N-glycosylated GLYCAM/OpenMM system

Use this guide when you have one or more residue-resolved, canonical
GLYCAM-named glycan PDB fragments and want PolyzyMD to attach them to Asn sites,
parameterize the product with native OpenMM GLYCAM, and run or export through the
PolyzyMD exact OpenMM/GROMACS paths.

```{important}
This is an explicit per-moiety route. Mark each glycan moiety with
`moiety.force_field: glycam06`. Moieties that omit `force_field` inherit the
top-level `force_field.small_molecule` setting and use the generic OpenFF route.
Unknown force-field labels fail loudly; PolyzyMD does not silently fall back from
GLYCAM to OpenFF.
```

```{warning}
In this workflow, **exact** means PolyzyMD owns exact force-field parameter
assignment and exact transfer of local OpenMM exception/exclusion semantics into
the PolyzyMD GROMACS export. It does not mean bitwise-identical full PME total
energies between OpenMM and GROMACS. PME mesh, interpolation order, tolerances,
cutoff modifiers, and related engine settings are engine hyperparameters that
users must configure deliberately and validate for their ensemble.
```

## Prerequisites

- A cleaned protein PDB containing the target Asn residue or residues with
  explicit hydrogens.
- One residue-resolved glycan PDB per attachment. The glycan must use canonical
  GLYCAM residue and atom names and must contain complete curated `CONECT`
  records for the glycan graph.
- The reducing-end anomeric `C1` must have explicit hydroxyl O/H atoms, either
  as an ordinary residue-local hydroxyl or as the supported separate-residue
  `ROH` cap with `O1`/`HO1` atom names.
- A PolyzyMD configuration that enables `conjugation` and defines one
  `n_glycosylation` attachment per Asn site using `moiety.input_path`.

The exact GLYCAM route supports multiple glycan attachments, for example Asn25
and Asn60 in the same protein. Each glycan attachment is scoped internally so
Pablo can parse repeated or branched residue names without collapsing identities;
PolyzyMD restores the canonical GLYCAM names after the Pablo-only matching step.

## 1. Find the glycan on RCSB PDB

1. Open the RCSB PDB structure page, for example
   <https://www.rcsb.org/structure/5FYJ>.
2. Scroll to the **Glycans** section.
3. Identify the glycan you want to model and note its GlyGen accession, such as
   `G80966KZ`.

## 2. Download or prepare a 3D residue-resolved glycan PDB

1. Open the GlyGen glycan page, for example
   <https://www.glygen.org/glycan/G80966KZ>.
2. Use the GlyGen accession to find a downloadable 3D PDB from GlyGen or another
   residue-resolved glycan source. RCSB identifies the glycan in the protein
   structure; it does not by itself provide the standalone reducing-end glycan
   PDB that this loader expects.
3. GlycoShape or other downloads may be usable after manual preparation: extract
   a single `MODEL`, retain residue labels, and verify explicit reducing-end
   C1-O-H chemistry.
4. Save it in your project, for example:

   ```text
   structures/G80966KZ_glycam.pdb
   ```

### File considerations

- The generic loader requires a single connected fragment graph and preserves
  residue labels.
- Files must include complete, curated `CONECT` records for the whole fragment.
  PolyzyMD rejects missing and detectably invalid explicit graphs, including
  unknown endpoints, self-bonds, disconnected or isolated atoms, invalid
  explicit-H degree, and obvious overvalence where applicable.
- PolyzyMD cannot prove every expected chemical bond is present and never
  repairs or infers omitted connectivity from coordinates. RDKit assigns bond
  orders only on the accepted exact connectivity. Inspect the output PDB and
  sidecar before using the structure downstream.
- Multi-model PDB files are not automatically reduced. If a GlycoShape or other
  download contains more than one `MODEL` record, extract the model you want into
  a single-model PDB before loading it.
- The input must contain a structurally valid reducing-end `C1` bonded to an
  explicit hydroxyl O/H group. The `ROH:O1`/`ROH:HO1` cap is accepted, but not
  required when an ordinary residue-local hydroxyl is present.

## 3. Configure force-field ownership

Set the ordinary protein and small-molecule force fields once. The GLYCAM route
is selected by `moiety.force_field: glycam06` on each glycan attachment, not by a
global glycan policy knob.

```yaml
force_field:
  protein: ff14sb_off_impropers_0.0.4.offxml
  small_molecule: openff-2.0.0.offxml
```

This route uses the OpenMM-vendored Amber14 XML stack for ff14SB protein,
GLYCAM_06j-1 glycan residues, GLYCAM NLN for modified Asn, and TIP3P water. The
native GLYCAM reference records its route invariants in the audit and exact
sidecar. In mixed GLYCAM/OpenFF overlays, the baseline/user simulation globals
such as PME method, cutoff, switching, and tolerance remain authoritative;
native-reference differences are diagnostics and do not change exact particle,
bonded, constraint, or exception transfer.
Disconnected precharged OpenFF components such as DMSO or ethanol, free polymer
molecules, and other disconnected organics can coexist with GLYCAM glycans when
they have complete provenance and assigned partial charges. Ions and water are
also included in the final solvated system. Covalently attached OpenFF polymer on
the same protein as a GLYCAM glycan is handled by the mixed overlay route, which
audits ownership and fails closed on unknown baseline parameter forces that touch
GLYCAM atoms.

Disconnected multi-residue Sage polymers are represented internally as one
OpenMM residue only for `SMIRNOFFTemplateGenerator` template matching. The
original monomer segmentation and component provenance remain in the audit
artifacts.

## 4. Configure one or more N-glycosylation attachments

Add one enabled attachment under `conjugation.attachments` for each Asn site. Use
the cleaned protein's Asn chain and residue number, set
`moiety.force_field: glycam06`, and point `moiety.input_path` at the canonical
GLYCAM-named CONECT PDB. The built-in N-glycosylation mechanism infers the
canonical Asn `ND2`--glycan reducing-end `C1` linkage and validates the hydroxyl
leaving atoms. Use explicit atom selectors only when you need to override that
canonical inference.

```yaml
conjugation:
  enabled: true
  ccd_pablo:
    enabled: true
    lookup_policy: auto_download
    use_canonical_atom_names: false
  attachments:
    - name: asn25_glycan
      site:
        chain_id: A
        residue_name: ASN
        residue_number: 25
        atom_name: ND2
      moiety:
        name: G42666_asn25
        force_field: glycam06
        input_path: structures/G42666_asn25_glycam.pdb
      mechanism:
        name: n_glycosylation
    - name: asn60_glycan
      site:
        chain_id: A
        residue_name: ASN
        residue_number: 60
        atom_name: ND2
      moiety:
        name: G80966_asn60
        force_field: glycam06
        input_path: structures/G80966_asn60_glycam.pdb
      mechanism:
        name: n_glycosylation
```

The built-in mechanism forms an Asn `ND2`--glycan reducing-end `C1` bond. It
removes one Asn `ND2` hydrogen and the validated glycan hydroxyl O/H leaving
atoms. PolyzyMD uses attachment-scoped internal residue aliases while Pablo parses
the full modified topology, then restores canonical glycan names and routes the
modified Asn residue to the GLYCAM `NLN` template during OpenMM system creation.

## 5. Use a tested YAML shape

The following configuration uses only implemented fields. Adjust file paths,
durations, and counts for your system.

```yaml
name: glycam_asn25_asn60
description: Strict GLYCAM N-glycosylated OpenMM build

enzyme:
  name: enzyme
  pdb_path: structures/enzyme_prepared.pdb

substrate: null

polymers:
  enabled: true
  generation_mode: cached
  type_prefix: FREEPOLY
  monomers:
    - label: A
      probability: 1.0
      name: EGPMA
      residue_name: EGP
  length: 4
  count: 2
  sdf_directory: structures/free_polymers

conjugation:
  enabled: true
  ccd_pablo:
    enabled: true
    lookup_policy: auto_download
    use_canonical_atom_names: false
  attachments:
    - name: asn25_glycan
      site:
        chain_id: A
        residue_name: ASN
        residue_number: 25
        atom_name: ND2
      moiety:
        name: G42666_asn25
        force_field: glycam06
        input_path: structures/G42666_asn25_glycam.pdb
      mechanism:
        name: n_glycosylation
    - name: asn60_glycan
      site:
        chain_id: A
        residue_name: ASN
        residue_number: 60
        atom_name: ND2
      moiety:
        name: G80966_asn60
        force_field: glycam06
        input_path: structures/G80966_asn60_glycam.pdb
      mechanism:
        name: n_glycosylation

solvent:
  primary:
    type: water
    model: tip3p
  co_solvents:
    - name: dmso
      mole_fraction: 0.10
      residue_name: DMS
  ions:
    neutralize: false
    nacl_concentration: 0.0
    kcl_concentration: 0.0
    mgcl2_concentration: 0.0
  box:
    padding: 1.2
    shape: rhombic_dodecahedron
    target_density: 1.0
    tolerance: 2.0

thermodynamics:
  temperature: 300.0
  pressure: 1.0

simulation_phases:
  equilibration_stages:
    - name: heat
      duration: 0.1
      samples: 10
      ensemble: NVT
      temperature: 300.0
      position_restraints:
        - group: protein_heavy
          force_constant: 4184.0
    - name: density
      duration: 0.1
      samples: 10
      ensemble: NPT
      temperature: 300.0
      barostat: MC
  production:
    ensemble: NPT
    duration: 1.0
    samples: 100
    report_interval: 5000
    time_step: 2.0
    checkpoint_interval: 600.0

output:
  projects_directory: outputs/projects
  scratch_directory: outputs/scratch
  naming_template: "{enzyme}_{duration}ns_{temperature}K_run{replicate}"
  trajectory_format: dcd

force_field:
  protein: ff14sb_off_impropers_0.0.4.offxml
  small_molecule: openff-2.0.0.offxml

engine: openmm
openmm:
  platform: CUDA
  precision: mixed
```

The final solvated system may contain the attached glycans, attached OpenFF
polymers, free OpenFF polymer chains, organic cosolvent, ions, and water. The
final periodic E2E acceptance test exercises this mixed component set.

## 6. Run the build

Use the normal PolyzyMD build entry point for your project. The config above
selects the native GLYCAM handoff internally; the returned conjugation result
contains an `ExactExportBundle` rather than a public vanilla Interchange.

```bash
pixi run -e build polyzymd build -c config.yaml
```

The OpenMM execution path uses the authoritative native OpenMM `System`,
`Topology`, and positions from that bundle. Production minimization should be
convergence-driven. Use uncapped OpenMM minimization (`maxIterations=0`) when
validating final exact-bundle stability rather than treating an arbitrary
iteration cap as convergence evidence.

## 7. Inspect provenance artifacts

The important exact-route artifacts are:

| Artifact | Typical path | What to check |
|----------|--------------|---------------|
| Solvated PDB | `.../solvated_conjugate.pdb` | Canonical glycan residue and atom names are restored after Pablo-only scoped aliasing. |
| Native GLYCAM audit | `.../native_openmm_glycam_audit.json` | `route: native_openmm_glycam`, `residue_templates` maps each modified Asn/ASX residue to `NLN`, `glycam_template_matches` lists matched GLYCAM residues, and `sage_template_generator.proof_no_glycan_entered_sage` is `true`. |
| Exact exception sidecar | `.../exact_openmm_exceptions.json` | Schema v2 sidecar with authoritative OpenMM atoms, bonds, nonbonded particles, exceptions, constraints, route invariants, and normalized GROMACS output identity/hashes. |
| Workflow JSON | `.../conjugated_polymer_system_workflow.json` | `artifact_paths.native_openmm_glycam_audit` and `artifact_paths.exact_openmm_exceptions` point to the audit and sidecar. |

For quick inspection:

```bash
python -m json.tool outputs/scratch/.../native_openmm_glycam_audit.json | less
python -m json.tool outputs/scratch/.../exact_openmm_exceptions.json | less
```

## 8. Export exact GROMACS files only through PolyzyMD

For `engine: gromacs`, or when using `GromacsExporter` with an exact bundle,
PolyzyMD writes a baseline topology through a private Interchange and then patches
it with the authoritative OpenMM exception sidecar. Before patching, the exporter
checks count, atom-order identity, topology-bond identity, and normalized GROMACS
hashes against the sidecar and fails closed on mismatch. The patched topology
sets `gen-pairs` to `no`, disables automatic bonded exclusions per molecule type,
adds explicit function-2 `[ pairs ]` rows for every nonzero OpenMM exception, and
adds exact `[ exclusions ]` rows. It then writes an audit such as
`<prefix>_exact_gromacs_audit.json`.

Do not call raw `Interchange.to_gromacs()` on the underlying private Interchange.
That raw export loses exact GLYCAM mixed scaled/unscaled explicit 1-4 exception
semantics.

After export, compare short OpenMM and GROMACS smoke trajectories or energy/force
checks under your chosen PME settings. Treat differences from PME grid/order or
modifier choices as engine configuration differences to manage, not as evidence
that PolyzyMD changed local bonded or exception parameters.

Exact GROMACS export currently does not generate glycan-specific position
restraints. If you need glycan restraints, use the exact OpenMM path with staged
post-overlay `position_restraints`, or provide and validate a separate GROMACS
restraint workflow outside the exact exporter.

The periodic final acceptance test is opt-in:

```bash
POLYZYMD_RUN_FINAL_E2E=1 pixi run -e build pytest tests/test_conjugation_final_e2e.py -m final_e2e -v
```

## Limitations

- The glycan input must be a canonical GLYCAM-named PDB with complete `CONECT`
  records. PolyzyMD rejects missing or detectably invalid explicit graphs and
  does not infer bonds from coordinates.
- Strict GLYCAM never silently falls back to OpenFF. Omit `moiety.force_field`
  only for moieties that should inherit `force_field.small_molecule`.
- Covalently attached Sage polymer on the same protein as a GLYCAM glycan is
  unsupported because cross-boundary bonded and exception provenance is not
  audited. Disconnected Sage components are supported when precharged and routed
  through `SMIRNOFFTemplateGenerator`.
- Pablo topology representability has degree and residue-template limits. If the
  full modified topology is not representable, PolyzyMD fails closed.
- Only PolyzyMD's exact OpenMM and exact GROMACS paths are guaranteed for this
  workflow. Arbitrary vanilla Interchange exporters are not guaranteed.
- Exact GROMACS export currently does not generate glycan-specific position
  restraints and rejects `component_info` position-restraint postprocessing. If
  glycan position restraints are required, use the exact OpenMM route with staged
  post-overlay `position_restraints`, or provide a separately validated GROMACS
  restraint workflow.

## Related reference

- {doc}`../reference/protein_modification_config`
- {doc}`../reference/conjugation_support_matrix`
- {doc}`../reference/glycam_exact_export`
- {doc}`../reference/openff_pdb_ingestion`
