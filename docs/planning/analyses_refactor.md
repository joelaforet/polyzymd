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
  - Status: in review

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
  - Status: not started

- [ ] Restrict hydrogen bond donors and acceptors to N and O, with S optional,
      and record the effective selections in provenance.
  - Branch: `analyses/hbond-selections`
  - Owner:
  - Status: not started

- [ ] Thread `ttest_method` and `posthoc_method` through every plugin
      `compare()`, add the Benjamini-Hochberg correction to distances, and
      define the correction family once so the ANOVA applies it consistently.
  - Branch: `analyses/test-threading`
  - Owner:
  - Status: not started

- [ ] Add `ci95_low`, `ci95_high`, `ci_method` and `unit` to every metric model,
      label every error bar and band with what it means and over how many
      replicates, and add an `uncertainty` block to every artifact.
  - Branch: `analyses/confidence-intervals`
  - Owner:
  - Status: not started

- [ ] Consult the per-segment status in `progress.json` before loading, and
      refuse a cached replicate whose recorded size or mtime differs from what
      is on disk.
  - Branch: `analyses/segment-cache-freshness`
  - Owner:
  - Status: not started

- [ ] Declare a periodic boundary policy on load, record it, and fail rather
      than warn when fragment mode or chain identity needs bonds that are
      absent. Drop alignment in distances and the catalytic triad, or disable
      the minimum image convention when the frames are aligned.
  - Branch: `analyses/pbc-alignment`
  - Owner:
  - Status: not started

## Phase 1, agent protocol (v1.3.0)

- [ ] Add the `polyzymd analyze` command, `protocols.analyze()`,
      `ProtocolReport` and `--format agent`, plus the `polyzymd-analyze` skill.
      Fix the seaborn import.
  - Branch: `analyses/agent-protocol`
  - Owner:
  - Status: not started

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
  - Status: not started

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
