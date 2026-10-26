# Evaluate PolyzyMD

This page is the fast path for a reviewer, or anyone who wants to judge
PolyzyMD before using it. It lists, in order, what to run and what to read,
and how long each step takes: about 45 minutes in all. Each step links the
page that holds the evidence.

Start with {ref}`scope-and-limitations`, which says in a few paragraphs what
PolyzyMD builds and runs and what it does not.

## 1. Run the quickstart (3 minutes)

Install the `build` environment ({doc}`index`) and follow
{doc}`quickstart`. It makes a project, builds Trp-cage in water with NaCl
(6,138 atoms), audits the build, runs a few picoseconds on the CPU and
reports a radius of gyration. On a laptop CPU, `polyzymd run` takes about 45
seconds.

## 2. Open the build record (2 minutes)

Open `build_provenance.md` in the replicate folder,
`runs/trpcage/water/trpcage_300K_run1/`. Every build writes it. It holds:

- the PolyzyMD and OpenMM versions, the {term}`build config hash` and the
  atom count;
- a table of components with their force field files and charge sources;
- the physics audit, one row per check and component, with PASS, WARN, FAIL
  or INFO and the numbers behind it (25 pass and 2 info for the quickstart);
- the methods to cite, with DOIs, for exactly what the build used;
- the seeds and box vectors.

`build_provenance.json` holds the same record, with the SHA-256 of every file
the engine reads.

## 3. Read how the physics audit works (10 minutes)

{doc}`../explanation/physics_audit` explains what the audit compares, how it
rebuilds each component on its own with OpenFF, what happens when a check
fails, how a deliberate departure is recorded, and what the audit does not
check. {doc}`../reference/physics_checks` lists every check.

## 4. See what is tested and where (5 minutes)

| Claim | Tests |
|---|---|
| Each audit check fails on a build broken in that way | `tests/audit/` (see {ref}`the table on the audit page <physics-audit-tests>`) |
| OpenMM and GROMACS exports hold the same parameters and settings | `tests/audit/test_components_and_engines.py`, `tests/audit/test_engine_parity.py` |
| GROMACS applies the staged position restraints PolyzyMD writes | `tests/exporters/test_gromacs_staged_restraints.py` |
| The quickstart runs end to end on both engines | `tests/test_quickstart_example.py` |
| Seeds depend only on the replicate and the phase | `tests/simulation/test_seeds.py` |
| Statistics: known answers for the statistical inefficiency, the intervals and the tests | `tests/analyses/scientific/`, `tests/analyses/shared/test_inferential_statistics.py` |
| Freezing a study or project | `tests/analyses/test_study_freeze.py`, `tests/analyses/test_project_freeze.py` |

The GitHub Actions workflows run them
([all runs](https://github.com/joelaforet/polyzymd/actions)):

- `.github/workflows/engine-parity.yml`: the audit and engine parity tests,
  with GROMACS installed, on every pull request.
- `.github/workflows/full-test.yml`: the whole test suite on Linux and macOS,
  for pull requests to `main` and `dev` and once a week.
- `.github/workflows/ci.yml`: lint, type check, package build, install test
  and a docs build with warnings as errors.

## 5. Compare the two engines (5 minutes)

{doc}`../reference/gromacs_openmm` gives, for each config setting, what
OpenMM and GROMACS run, and marks each row 1:1, approximate or differs.
{doc}`../reference/validation` gives what the CI compares between the
engines, and a single-point energy of the quickstart system on both engines,
term by term, with the script to repeat it.

## 6. Check the statistics (5 minutes)

{doc}`../explanation/analysis_statistics_best_practices` explains how
PolyzyMD reports uncertainty: one value per replicate, a Student-t interval,
Welch's test between conditions and Benjamini-Hochberg across a family of
tests. It also lists what PolyzyMD does not estimate. To see an interval, run
two more replicates of the quickstart (about 45 seconds each on a CPU) and
analyze again:

```bash
polyzymd run -c trpcage/conditions/water/config.yaml -r 2-3
polyzymd analyze --project .
```

## 7. Read the verification records (5 minutes)

Each record compares one analysis with an independent calculation on test
systems, and names its tests:
{doc}`../reference/analysis_rmsf_verification`,
{doc}`../reference/analysis_sasa_verification`,
{doc}`../reference/analysis_contacts_verification`,
{doc}`../reference/analysis_native_contacts_verification` and
{doc}`../reference/analysis_hydrogen_bonds_verification`.

## 8. Check reproducibility (5 minutes)

- **Seeds.** The replicate number seeds the starting structure, the
  velocities and the thermostat and barostat noise of each stage. Runs are
  not bitwise reproducible. See {doc}`../reference/gromacs_openmm` (Seeds).
- **Freeze and reproduce.** `polyzymd project freeze` writes a deposit with a
  SHA-256 manifest, `CITATION.cff` and the engine inputs, and tags the commit.
  A reader checks the download with `polyzymd study locate --verify` and
  recomputes the analyses with `polyzymd analyze --recompute`. See
  {doc}`../how_to/study_freeze`. Freezing commits and tags, so try it on a
  copy of the quickstart project.

## 9. References and citation (2 minutes)

{doc}`../reference/references` lists the methods PolyzyMD implements and the
software it uses. Each build lists its own methods to cite in
`build_provenance.md`. To cite PolyzyMD, use `CITATION.cff` in the
repository.

## 10. Scope and limitations

Read {ref}`scope-and-limitations` again with what you have seen: the
engines, force fields, charge models, ensembles and statistics that PolyzyMD
supports, and what is out of scope.
