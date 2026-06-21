# Validate Conjugate Build Artifacts

Use this guide after building an opt-in conjugated system to decide whether the
construction artifacts are internally consistent enough to inspect, export, or
carry forward into simulation setup.

Validation is a reliability check, not a scientific guarantee. A passing report
means PolyzyMD found the expected product connectivity, charge evidence,
parameterization evidence, and smoke-test evidence for the current supported
workflow. It does not prove that an arbitrary chemistry, force field choice, or
production protocol is physically correct.

## Before you start

You need a `config.yaml` with `conjugation.enabled: true` and at least one
enabled attachment. The validated executable path is the current config-driven
NHS-lysine polymer vertical slice:

- unmodified protein input;
- NHS-activated polymer recipe attached to lysine;
- Pablo/OpenFF Interchange product-state ingestion;
- ff14SB protein parameters plus polymer templates;
- local NAGL patch charge bridge around the linkage;
- local charge reconciliation;
- restrained OpenMM smoke/minimization evidence;
- final solvation.

The N-glycosylation path is wired for mechanism tests and exploratory workflow
development, but it is not the same validated vertical slice as the NHS-Lys
polymer workflow. See {doc}`../reference/conjugation_support_matrix` for the
current support levels.

## Build the conjugated system

Run the normal build in a full PolyzyMD build environment:

```bash
pixi run -e build polyzymd build -c config.yaml
```

Do not run full conjugation builds in lean simulation-only environments such as
`sim-cuda-12-4`. Those environments are intended for already prepared runtime
inputs. See {doc}`hpc_cuda_conjugation` for the two-stage HPC pattern and its
current handoff caveats.

## Find the validation report

The build writes a canonical validation report named:

```text
conjugate_validation_report.json
```

For public Python workflows, the report path is also included in
`ConjugationResult.artifact_paths` under `conjugate_validation_report` when the
construction object is available.

```python
report_path = result.artifact_paths.get("conjugate_validation_report")
print(report_path)
```

The report lives with the conjugate construction artifacts. In current
config-driven builds, related evidence sidecars include files such as
`product_state_charge_bridge.json`,
`product_state_charge_bridge_local_reconciliation.json`, `vacuum_smoke.json`,
`restrained_smoke_diagnostics.json`, and `pre_smoke_geometry.json` when those
checks ran.

## Read the top-level status

Open the report and inspect the top-level `status` first:

```python
import json
from pathlib import Path

# After result = build_conjugate_from_config(...) or build_conjugate(...)
report_path = result.artifact_paths["conjugate_validation_report"]
report = json.loads(Path(report_path).read_text())
print(report["status"])
```

If you are reading artifacts outside the Python result object, use the report
inside the construction artifact directory, for example:

```python
report_path = Path("artifacts/lysine-polymer-conjugate/conjugate-construction/conjugate_validation_report.json")
report = json.loads(report_path.read_text())
```

Statuses are:

| Status | Meaning |
|--------|---------|
| `pass` | All available validation categories passed. |
| `warn` | No category failed, but at least one category found a conservative warning. Inspect before continuing. |
| `fail` | At least one required validation category failed. Do not treat the conjugate as ready. |
| `skipped` | Metadata or evidence was unavailable for all categories. This is not a validation pass. |

## Inspect each validation category

The report is organized into named audit sections. For each section, inspect
`status`, `checks`, `message`, and `evidence`.

| Report section | What it checks | Evidence to inspect |
|----------------|----------------|---------------------|
| `bond_graph` | Expected product linkage bonds are present in PDB `CONECT` records. | `expected_bonds`, `observed_bonds`, `missing_bonds`. |
| `atom_presence` | Required link atoms remain present and declared leaving atoms are absent. | `present_atoms`, `missing_atoms`, `lingering_leaving_atoms`. |
| `valence_sanity` | Conservative bond-count sanity near linkage atoms. | `bond_counts` and any high-count warning evidence. |
| `charge_audit` | Product-state charge bridge, total/formal charge reconciliation, and local reconciliation sidecar status. | `bridge_report_path`, `reconciliation_report_path`, `total_charge_e`, `formal_charge_e`, correction fields. |
| `parameter_coverage` | Final OpenFF Interchange can be converted to an OpenMM system and, when known, has the expected particle count. | `expected_particle_count`, `observed_particle_count`. |
| `linkage_geometry` | Linkage distances and severe nonbonded close contacts in the product coordinates. | `linkage_distances_angstrom`, `close_contact_count`. |
| `openmm_smoke` | OpenMM smoke and pre-smoke JSON evidence report success and finite diagnostics. | `smoke_json_path`, `diagnostics_json_path`, `pre_smoke_geometry_json_path`. |

For quick triage, print any non-passing sections:

```python
for section, payload in report.items():
    if isinstance(payload, dict) and payload.get("status") not in {"pass", None}:
        print(section, payload.get("status"))
        for check in payload.get("checks", []):
            print("  -", check.get("name"), check.get("message"))
```

## Decide what to do next

- Continue only after you understand every `warn` and have no `fail` statuses.
- Treat `skipped` as missing evidence, not as success, unless you intentionally
  disabled or bypassed that part of the workflow.
- If charge or parameter coverage fails, inspect the charge bridge and Pablo or
  Interchange diagnostics before rerunning expensive simulations.
- If geometry or smoke fails, inspect the crosslinked and minimized PDBs in a
  molecular viewer before changing simulation settings.

## Known limits

- Validation is necessary but not sufficient. It catches many construction
  failures, but it does not certify arbitrary chemistry or force-field accuracy.
- Current validated reliability work centers on the NHS-Lys polymer path. Other
  mechanisms may be wired or exploratory rather than production-ready.
- The exact Interchange/runtime handoff for CUDA 12.4 simulation-only
  environments is planned and implementation-dependent. Use only commands that
  appear in `polyzymd --help` or {doc}`../reference/cli_reference`.
- Platform fail-fast behavior for incompatible OpenMM/CUDA environments is a
  follow-up item; run explicit environment checks before submitting production
  jobs.
