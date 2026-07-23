# Process a Multi-Segment Protein with File-Based Attachments

This tutorial is a workflow template for building one prepared, multi-segment protein with several
file-based covalent attachments. It uses an attachment CSV so the input files
can stay in a separate folder, then stops at each public build scope for
inspection.

You will create four kinds of preparation artifact:

1. an assembled coordinate structure;
2. a parameterized, solvent-free solute;
3. optionally, a topology-only GROMACS handoff for that solute; and
4. a complete solvated system.

No command in this tutorial starts molecular dynamics. The optional GROMACS
command writes handoff files; it does not invoke GROMACS.

```{important}
The built-in executable file-based route in this tutorial is strict
N-glycosylation. The fragment filenames and number of attachments are your
choice, but each fragment must be a chemically suitable, residue-resolved
N-glycan PDB with canonical GLYCAM names and a complete, curated `CONECT`
graph. A folder of arbitrary PDB chemistries is not accepted automatically.
See {doc}`../reference/protein_modification_config` for the current mechanism
and input contract.
```

## Prerequisites

You need:

- PolyzyMD installed through the repository's `build` pixi environment;
- one prepared protein PDB containing all intended peptide segments and stable
  chain, residue, and atom labels;
- `TER` records at every peptide-segment boundary;
- explicit source `CONECT` records for covalent links that must survive, such
  as cross-segment disulfides;
- one complete-`CONECT` N-glycan PDB per attachment; and
- enough memory and wall time for parameterization and solvation.

The protein must already have the intended heavy atoms, protonation decisions,
and segment boundaries. PolyzyMD does not decide the biological assembly for
you. The default conjugation preparation removes source hydrogens and regenerates
them at pH 7, so do not rely on source hydrogen names or serials as selectors.
The fragment PDBs must describe their complete bond graphs. PolyzyMD does not
infer or repair omitted fragment bonds from coordinates.

Use {doc}`prepare_pdb_for_openff` if the protein still needs preparation. For
the detailed protein and fragment requirements, read
{doc}`../reference/openff_pdb_ingestion` and
{doc}`../how_to/troubleshoot_openff_pdb_ingestion`.

## Step 1: Arrange the input files

Create this project layout:

```text
multisegment-build/
|- config.yaml
|- structures/
|  `- prepared_multisegment.pdb
|- fragments/
|  |- attachment_01.pdb
|  `- attachment_02.pdb
`- tables/
   `- attachments.csv
```

Change into the project before running the remaining commands:

```bash
cd multisegment-build
```

This matters because the relative `output.projects_directory` below is resolved
from the command's working directory, not from the location of `config.yaml`.

The fragment filenames are only labels. PolyzyMD does not scan the directory;
every file to process must have one row in `attachments.csv`.

Before building, count the records that carry segment and connectivity
information:

```bash
grep -c '^TER' structures/prepared_multisegment.pdb
grep -c '^CONECT' structures/prepared_multisegment.pdb
grep -c '^CONECT' fragments/*.pdb
```

A nonzero count is not proof that the chemistry is correct. Open the protein
and every fragment in a molecular viewer. Confirm the intended segment breaks,
cross-segment bonds, reducing-end atoms, explicit hydrogens, and complete
fragment connectivity.

## Step 2: Describe attachments with source labels

Suppose the prepared protein has target Asn residues on source chains `A` and
`B`. Write `tables/attachments.csv` as:

```text
name,site.chain_id,site.residue_name,site.residue_number,site.atom_name,moiety.name,moiety.force_field,moiety.input_path,mechanism.name
site_01,A,ASN,25,ND2,attachment_01,glycam06,../fragments/attachment_01.pdb,n_glycosylation
site_02,B,ASN,60,ND2,attachment_02,glycam06,../fragments/attachment_02.pdb,n_glycosylation
```

Use the chain ID, residue name, residue number, and atom name from the prepared
input PDB. Do not use atom serial or atom index selectors for this workflow.
PolyzyMD's default protein-preparation step canonicalizes output to chain `A`
and does not expose an original-to-prepared identity ledger. Chain and residue
labels must therefore survive preparation well enough for the configured site
to resolve. The structure-scope build is the required check: verify every
resolved site and preserved `TER` boundary in its output instead of assuming
that an input label survived.

The CSV path is relative to `config.yaml`. Each `moiety.input_path` is relative
to the CSV file, which is why these rows use `../fragments/`.

CSV rows accept scalar fields only. Use inline `conjugation.attachments` in YAML
if a mechanism needs list-valued leaving atoms or other nested chemistry. The
complete header list and path rules are in
{doc}`../reference/protein_modification_config`.

## Step 3: Write a minimal valid configuration

Create `config.yaml`:

```yaml
name: multisegment_file_build
engine: openmm

build:
  scope: system

enzyme:
  name: prepared_protein
  pdb_path: structures/prepared_multisegment.pdb

conjugation:
  enabled: true
  attachments_file: tables/attachments.csv
  placement:
    timeout_seconds: 900

thermodynamics:
  temperature: 300.0

simulation_phases:
  equilibration_stages:
    - name: inspect_only
      duration: 0.01
      samples: 2
      ensemble: NVT
      temperature: 300.0
  production:
    ensemble: NPT
    duration: 0.01
    samples: 2
    checkpoint_interval: 60.0

output:
  projects_directory: outputs
  naming_template: "replicate_{replicate}"
```

The phase fields make this a schema-complete `SimulationConfig`; this tutorial
does not execute those phases. Default solvent, ion, box, and force-field
settings fill the omitted sections. Review and replace those defaults before
using the final system in a scientific simulation.

Validate the YAML and attachment table:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd build -c config.yaml --scope structure --dry-run
```

The dry run should report `Scope: structure`, the conjugation workflow, and a
checkpoint before parameterization. Fix referenced-file warnings before the
real build. Validation and dry run check the schema and referenced paths only.
They do not parse protein chemistry, prove a fragment's complete `CONECT` graph,
find its reducing-end motif, or resolve an attachment selector. Those checks
first occur in the real structure build.

## Step 4: Build and inspect the coordinate structure

Build only the covalently assembled coordinates:

```bash
pixi run -e build polyzymd build -c config.yaml --scope structure
```

The first checkpoint is:

```text
outputs/replicate_1/
|- conjugate-construction/
|  |- assembled_crosslinked.pdb
|  `- packmol_modifier_placement/
|     `- packmol_run_status.json
`- conjugated_polymer_system_workflow.json
```

Retry directories may instead be named
`packmol_modifier_placement_attempt_XX`. Exact diagnostic paths are recorded in
the workflow JSON.

Inspect `assembled_crosslinked.pdb` in a molecular viewer. Check all of the
following before parameterization:

- every intended attachment is present exactly once;
- every attachment is bonded to the intended Asn;
- each original protein segment boundary still has a `TER` record;
- no peptide bond was invented across a `TER` boundary;
- intended source `CONECT` links that survived graph surgery remain present;
- the final protein--fragment and fragment--fragment geometry has no severe
  nonbonded clashes; and
- modified Asn is represented as the expected `NLN` product with one `HD21`.

For a text audit, inspect the workflow and placement status:

```bash
pixi run -e build python -m json.tool \
  outputs/replicate_1/conjugated_polymer_system_workflow.json | less
pixi run -e build python -m json.tool \
  outputs/replicate_1/conjugate-construction/packmol_modifier_placement/packmol_run_status.json | less
```

The placement status must say `accepted_after_validation: true`. The validator
checks atom count and order, finite coordinates, linkage length, final
`CONECT`, and graph-aware clashes. The mere existence of `packmol.pdb` is not a
success condition.

## Step 5: Build and audit the isolated solute

Parameterize only the protein and its covalently attached fragments:

```bash
pixi run -e build polyzymd build -c config.yaml --scope solute
```

This adds the isolated OpenMM handoff:

```text
outputs/replicate_1/solute/
|- solute.pdb
|- system.xml
|- openmm_build_audit.json
|- native_openmm_glycam_audit.json
`- exact_openmm_exceptions.json
```

Inspect `solute.pdb` visually again. Confirm atom order and particle counts in
`openmm_build_audit.json`, then check that the native GLYCAM audit records the
intended route and that no glycan entered Sage:

```bash
pixi run -e build python -m json.tool \
  outputs/replicate_1/solute/openmm_build_audit.json | less
pixi run -e build python -m json.tool \
  outputs/replicate_1/solute/native_openmm_glycam_audit.json | less
```

The exact exception sidecar is authoritative parameter evidence. Do not edit it
or detach it from the coordinate/topology artifacts whose identities it records.

Solute scope intentionally excludes substrate, free polymer, water, ions,
barostats, configured stage restraints, and a runnable dynamics protocol. A
charged exact solute may retain a periodic PME representation to preserve
parameters. That does not make it a physically complete periodic system.

## Step 6: Optionally create a topology-only GROMACS handoff

If another workflow needs the isolated conjugate in GROMACS format, run:

```bash
pixi run -e build polyzymd build \
  -c config.yaml \
  --scope solute \
  --format gromacs
```

The handoff appears under:

```text
outputs/replicate_1/solute/gromacs/
|- *.gro
|- *.top
|- *.itp
|- *exact_gromacs_audit.json
`- exact_openmm_exceptions.json
```

Some topologies do not need split `.itp` files. Use the emitted topology and
audit as the source of truth rather than requiring a fixed count.

This command does **not** run `gmx`, energy minimization, equilibration, or
production. It deliberately emits no `.mdp`, run script, or position-restraint
artifact for solute scope. The receiving workflow must provide a physical box,
solvent, ions, electrostatics settings, restraints if needed, and an independently
validated run protocol. See {doc}`../how_to/gromacs_export` for the distinction
between a build-only handoff and an executed GROMACS workflow.

## Step 7: Build the complete solvated system

After the structure and solute audits pass, build the configured system scope:

```bash
pixi run -e build polyzymd build -c config.yaml --scope system
```

The final build adds route-specific evidence and the standard OpenMM inputs,
including:

```text
outputs/replicate_1/
|- conjugate-construction/
|  |- assembled_crosslinked.pdb
|  |- conjugate_relaxed.pdb
|  `- conjugate_validation_report.json
|- conjugated_polymer_system_workflow.json
|- native_openmm_glycam_audit.json
|- exact_openmm_exceptions.json
|- solvated_conjugate_free_polymers.pdb
|- solvated_system.pdb
`- system.xml
```

The route-specific solvated PDB and `solvated_system.pdb` alias should describe
the same final coordinates. Review the validation report, workflow JSON, native
GLYCAM audit, exact exception sidecar, solvent and ion composition, periodic box,
and OpenMM particle count. Building these files is preparation evidence, not a
claim that the molecular model is stable or scientifically validated.

The exact artifact tree can include more diagnostics than shown here. The
canonical scope table in {ref}`build-scope-artifacts` explains which files each
public endpoint promises.

## Understand timeout and Packmol exit 173

`conjugation.placement.timeout_seconds` is a per-process limit. Its default is
900 seconds. Set a shorter positive value while diagnosing expensive placement,
or set it to `null` only when you intentionally accept an unbounded wait.

A timeout means PolyzyMD terminated that Packmol process. Exit code 173 means
Packmol returned an imperfect best candidate. Neither result means the build
succeeded. PolyzyMD may continue only if the candidate passes the mandatory
validation recorded in `packmol_run_status.json`.

Single-fragment coordinate placement can make up to 20 sequential attempts, so
total wall time can exceed one timeout interval. A multi-fragment request uses
one joint attempt. If `accepted_after_validation` is false, treat the candidate
as rejected and fix the geometry or connectivity. The troubleshooting details
are in {doc}`../how_to/troubleshoot_openff_pdb_ingestion`.

## Know the staged-build limitation

The scopes are inspection endpoints, not a resumable pipeline. The current CLI
cannot start a solute or system build from a previously written
`assembled_crosslinked.pdb`. Each `polyzymd build` invocation reads the original
protein, attachment CSV, and fragment files and runs the required preceding
work again.

Keep all source files unchanged between scope checks. If you change a source
PDB, CSV row, YAML field, or fragment graph, repeat the earlier visual and audit
checks. Do not copy an assembled file over the source protein to simulate a
resume operation.

## What you should have now

You have learned how to:

- preserve source segment identity while selecting attachment sites;
- enumerate a folder of complete-`CONECT` fragments with a schema-checked CSV;
- inspect the structure, solute, optional GROMACS, and system endpoints; and
- distinguish generated preparation evidence from an executed or validated
  simulation.

For exact GLYCAM semantics and limitations, continue with
{doc}`../how_to/glygen_glycan_conjugation` and
{doc}`../reference/glycam_exact_export`. For command options, see
{doc}`../reference/cli_reference`.
