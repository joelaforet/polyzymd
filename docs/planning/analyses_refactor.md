# Analyses refactor checklist

This is the working checklist for the refactor of `src/polyzymd/analyses/`. It
comes from section 9 of the audit in
[analyses_audit_2026-09-11.md](analyses_audit_2026-09-11.md), which holds the
evidence for every item.

This file is not part of the Sphinx documentation. It is a plan, and it changes
as the work lands.

## How to use it

One item is one branch, one session and one pull request. Branch from
`analyses_refactor` using the name given under the item, and open the pull
request against `analyses_refactor`. Never commit on `analyses_refactor`
itself.

Fill in the owner when you pick an item up. Status is one of `not started`,
`in progress`, `in review`, or `merged`. Check the box when the pull request
merges, and do it in the same session that finishes the work so a crash costs
one item rather than the plan.

## Phase 0, correctness (v1.3.0)

- [ ] Make the repository legible to a fresh session. Add `CLAUDE.md` importing
      `AGENTS.md`, rewrite `.opencode/instructions/analysis-module.md` from the
      real module tree, add this checklist, add the `livecoms-check` skill, add
      the commit gate hook, and add `CITATION.cff` and the references page.
  - Branch: `analyses/infra`
  - Owner:
  - Status: in review (PR 105)

- [ ] Make `estimate_correlation_time` a thin wrapper around
      `statistical_inefficiency()`, remove the `max(tau, dt)` floor, and re-pin
      the tests to g near 1 on white noise and g near (1+phi)/(1-phi) on AR(1).
      Stop subsampling RMSF frames, keep tau as a diagnostic, and put
      uncertainty on RMSF through an across-replicate SEM. Fix the contacts
      aggregator docstring, which promises an autocorrelation correction it does
      not compute, and delete its Flyvbjerg citation until block averaging
      exists.
  - Branch: `analyses/correlation-time`
  - Owner:
  - Status: in review (PR 103)

- [ ] Restrict hydrogen bond donors and acceptors to N and O, with S optional,
      and record the effective selections in provenance.
  - Branch: `analyses/hbond-selections`
  - Owner:
  - Status: in review (PR 102)

- [ ] Thread `ttest_method` and `posthoc_method` through every plugin
      `compare()`, add the Benjamini-Hochberg correction to distances, and
      define the correction family once so the ANOVA applies it consistently.
  - Branch: `analyses/test-threading`
  - Owner:
  - Status: in review (PR 104)

- [ ] Add `ci95_low`, `ci95_high`, `ci_method` and `unit` to every metric model,
      label every error bar and band with what it means and over how many
      replicates, and add an `uncertainty` block to every artifact.
  - Branch: `analyses/confidence-intervals`
  - Owner:
  - Status: in review (PR 109, data layer; PR 110 `analyses/confidence-interval-plots` stacked on it)

- [ ] Consult the per-segment status in `progress.json` before loading, and
      refuse a cached replicate whose recorded size or mtime differs from what
      is on disk.
  - Branch: `analyses/segment-cache-freshness`
  - Owner:
  - Status: in review (PR 106; PR 107 `analyses/cache-freshness` and PR 108 `analyses/pbc-alignment` stacked on it)

- [ ] Declare a periodic boundary policy on load, record it, and fail rather
      than warn when fragment mode or chain identity needs bonds that are
      absent. Drop alignment in distances and the catalytic triad, or disable
      the minimum image convention when the frames are aligned.
  - Branch: `analyses/pbc-alignment`
  - Owner:
  - Status: in review (PR 108, stacked on PR 106)

## Phase 1, agent protocol (v1.3.0)

- [ ] Add the `polyzymd analyze` command, `protocols.analyze()`,
      `ProtocolReport` and `--format agent`, plus the `polyzymd-analyze` skill.
      Fix the seaborn import.
  - Branch: `analyses/agent-protocol`
  - Owner:
  - Status: in review (PR 111, stacked on PR 109)

- [ ] Add `load(config, replicate, window=...)` and
      `load_files(topology, trajectories, dt_ps=...)` returning a
      `LoadedTrajectory` with a `LoadProvenance`, and write down in one place
      which loader is canonical.
  - Branch: `analyses/loader-api`
  - Owner:
  - Status: not started

## Phase 2, trust (v1.3.0)

- [ ] Write the known-answer tests under `tests/analyses/scientific/`. White
      noise and AR(1) for both correlation estimators, SEM calibration over
      synthetic replicates, Kabsch on a rotated cloud, single-sphere SASA
      against 4 pi (r+p) squared, DSSP on an ideal helix, Rg on a cube, contacts
      on a two-atom toy, hydrogen bonds on a toy with one N-H to O pair and one
      C-H to O pair, and Welch against Student on unequal-variance fixtures.
      Commit one tiny real trajectory, at most 50 frames and 500 atoms, that
      exercises equilibration and stride end to end and cross-checks RMSF
      against `MDAnalysis.analysis.rms.RMSF`.
  - Branch: `analyses/scientific-tests`
  - Owner:
  - Status: not started

- [ ] Add a NumPy-style `References` section and a one-line method statement to
      every plugin `__init__.py`, print the method statement under the results
      table in the default formatter, and fix the stale reference pages listed
      in section 8 of the audit. `CITATION.cff` and
      `docs/source/explanation/references.md` landed with `analyses/infra`.
  - Branch: `analyses/references-docs`
  - Owner:
  - Status: not started

- [ ] Rewrite the contributor guide so it matches the scaffold it documents,
      covers `extract_metrics()`, and asks a contributor to cite the method,
      write a known-answer test, and state what the error bar means.
  - Branch: `analyses/contributor-guide`
  - Owner:
  - Status: not started

## Phases 3 to 5, framework collapse (v1.4.0)

Run these in order. Each one depends on the `Observable` contract landing
first, so do not fan them out.

- [ ] Introduce `Observable(name, unit, kind, values, index)` with `kind` one of
      `mean_of_timeseries`, `fluctuation`, `fraction`, `distribution` or
      `profile`, reduce a plugin to `Settings` plus
      `compute(universe, frames, settings) -> Sequence[Observable]`, and write
      the contract tests.
  - Branch: `analyses/observable-contract`
  - Owner:
  - Status: in review (PR 112; `rg2` prototype, `--style contract` scaffold, `polyzymd-extend` skill)

- [ ] Give the framework one persistence path with a framework-written identity
      block holding the polyzymd version, the plugin source hash, the settings
      fingerprint, the config hash and the input file identity, and key
      replicate cache reuse on that block.
  - Branch: `analyses/framework-persistence`
  - Owner:
  - Status: not started

- [ ] Give the framework one aggregation path per `kind`, covering the
      N_eff-corrected mean and SEM, fluctuation, binomial fraction, kernel
      density estimate and per-index profile.
  - Branch: `analyses/framework-aggregation`
  - Owner:
  - Status: not started

- [ ] Give the framework one comparison path, taking Welch or Student and
      Benjamini-Hochberg or Tukey from a single defaults object, and one
      generic plotting and formatting path per `kind` with an optional
      `extra_plots` hook.
  - Branch: `analyses/framework-comparison`
  - Owner:
  - Status: not started

- [ ] Collapse the `MDA*Context` classes into the four framework contexts and
      delete the `ConditionSummary` and `ComparisonResult` family in favour of
      the artifact envelope.
  - Branch: `analyses/context-collapse`
  - Owner:
  - Status: not started

- [ ] Remove the `__module__` rewriting and the removed-hook police, give
      `Analysis` real abstract methods or make it a Protocol, and replace
      `SimulationConfig` in contexts with a `TrajectorySource` protocol.
  - Branch: `analyses/plugin-protocol`
  - Owner:
  - Status: not started

- [ ] Add the `polyzymd.analyses` entry-point group so out-of-tree plugins load
      without living in the package.
  - Branch: `analyses/entry-points`
  - Owner:
  - Status: not started

- [ ] Port rmsf, rg, rmsd and sasa to `Observable`, rmsf first because it is
      smallest and already declares a statistical policy.
  - Branch: `analyses/port-scalar-plugins`
  - Owner:
  - Status: not started

- [ ] Port the plotters to the generic per-kind plotting path.
  - Branch: `analyses/port-plotters`
  - Owner:
  - Status: not started

- [ ] Port contacts, hydrogen bonds and the catalytic triad, adding an optional
      `Aggregator` protocol for residence times.
  - Branch: `analyses/port-event-plugins`
  - Owner:
  - Status: not started

## What stays untouched

The audit asks that these keep their current shape through the collapse:
`mda/artifacts.py`, `mda/store.py`, `mda/frame_selection.py`, `mda/job.py`,
`mda/universe.py`, `shared/inferential_statistics.py`,
`shared/autocorrelation.py` after the estimator fix, the internals of
`shared/loader.py`, `discovery.py` and `exceptions.py`.

## Source-line ledger

Joe's standing requirement is that this refactor reduces lines and complexity. Every pull request reports its net change in source lines (files under `src/`, tests and docs excluded) and the reviewer blocks unexplained growth. Correctness fixes carry tests and a few new fields, so Phase 0 grows the source; the reduction comes from the Observable contract and the plugin ports that follow it, which delete the per-plugin aggregation, comparison, formatting and plotting stacks.

| PR | Branch | Source added | Source removed | Net |
|---|---|---|---|---|
| 105 | `analyses/infra` | 0 | 0 | 0 |
| 102 | `analyses/hbond-selections` | 304 | 35 | +269 |
| 103 | `analyses/correlation-time` | 263 | 288 | -25 |
| 104 | `analyses/test-threading` | 857 | 352 | +505 |
| 106 | `analyses/segment-cache-freshness` | 288 | 36 | +252 |
| 107 | `analyses/cache-freshness` | 594 | 69 | +525 |
| 108 | `analyses/pbc-alignment` | 570 | 135 | +435 |
| 109 | `analyses/confidence-intervals` | 913 | 408 | +505 |
| 110 | `analyses/confidence-interval-plots` | 660 | 87 | +573 |
| 111 | `analyses/agent-protocol` | 1442 | 22 | +1420 |
| 112 | `analyses/observable-contract` | 1357 | 6 | +1351 (rg reimplemented in 84 lines against 3,881; the eight ports that follow delete an estimated 10,000 to 12,000) |

Figures are from `git diff --numstat <base>..<tip> -- 'src/**/*.py'` at the time each pull request was last reviewed; update the row when a branch changes.

## Validation on real data

The LipA 363 K campaign (`/projects/jola3134/Enzyme_Immobilization/polyzymd_sims_config_and_run_files/LipA_363K_REDO`, six conditions, five replicates) was used read-only to check the Phase 0 branches against real artifacts and topologies. In the 50:50 replicate 1 hydrogen bond artifact only 0.8 percent of 16.1 million recorded events had nitrogen or oxygen at both ends, and the reported protein-polymer count of 132.9 per frame falls to 12.4 when carbon and self-paired hydrogen events are removed. Every campaign RMSF profile was computed from 4 of 2,000 production frames. Polymer fragment-mode Rg used 38 real fragments with no fallback. Real runs keep one `progress.json` at the run root and interrupted segments are the normal restart case. The full notes are kept with the audit outside the repository.

## Port phase (started 12 September 2026)

Joe asked for the plugin ports to start before the eleven pull requests above are reviewed. The branch `analyses_ported` merges PRs 102 to 112 as they stand and is the base for every port; conflicts with the eventual merges into `analyses_refactor` are dealt with afterwards. Each port is one branch `analyses/port-<plugin>` and one pull request against `analyses_ported`. A port reimplements the plugin on the Observable contract following `rg_contract` and the `--style contract` scaffold, freezes reference values from the old implementation on real LipA data (two runs copied to local storage; frames 500 to 600), proves parity to the frozen values, deletes the old package outright, and then deletes whatever that leaves unreferenced. The prerequisite `analyses/contract-plots` supplies figures keyed on kind so the ports can delete their `_plotters.py` modules.

Order: wave 1 is contract plots plus rmsf, secondary structure, rmsd and sasa; wave 2 is rg (promoting `rg2`), distances and catalytic triad, which share the pair-distance machinery; wave 3 is hydrogen bonds and contacts, the two with event-level outputs; wave 4 removes the legacy framework (`_framework/compare.py`, most of `comparison_models.py`, `mda/plugin.py`, `mda/comparison.py`, `_framework/contract.py`, the multi-run comparison and formatting helpers, the bespoke parts of `stats.py`, the `MDA*Context` classes, the `__module__` rewriting in `base.py`, the non-contract scaffold templates) and collapses `protocols.py` to one result shape. Target after wave 4: the module under 20,000 lines, from 51,982 at the start.

### Port phase pull requests (against `analyses_ported`)

| PR | Branch | Plugin or item | Source net | Parity on real data |
|---|---|---|---|---|
| 113 | `analyses/contract-plots` | figures keyed on kind | +460 | not applicable |
| 114 | `analyses/fix-alignment-reference` | `AlignTraj` reference frame bug, shared versions in identity | -62 | RMSF centroid and frame modes now differ (0.8437 vs 0.8434 Å) |
| 115 | `analyses/port-secondary-structure` | secondary_structure, 1,790 to 210 lines | -1,551 | exact |
| 116 | `analyses/port-rmsd` | rmsd, 3,131 to 304 lines; convergence module deleted | -3,025 | exact (average), float32 noise after dropping the alignment pass (frame, centroid, external) |
| 117 | `analyses/port-sasa` | sasa, 4,538 to 338 lines | -4,157 | exact, four contexts |
| 118 | `analyses/port-rmsf` | rmsf, 2,195 to 273 lines | -1,797 | exact, three modes |
| 119 | `analyses/port-distances` | distances and catalytic_triad, 5,892 to 326 lines | -6,025 | exact, five pairs |
| 120 | `analyses/port-rg` | rg, 3,987 to 288 lines | -3,555 | exact, protein and 38 fragments |
| pending | `analyses/port-hydrogen-bonds` | hydrogen_bonds, 4,097 to 430 lines | -3,596 | exact, seven partitions; corrected protein-polymer count 7.07 per frame against 132.9 |
| pending | `analyses/port-contacts` | contacts, 6,203 to 420 lines; shared/groupings and shared/selectors deleted | -7,378 | exact |

Sum of the ports so far: about 30,700 source lines removed against about 1,900 added, before the legacy framework removal.
