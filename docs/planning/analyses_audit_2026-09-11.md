# PolyzyMD analyses module audit for v1.3.0

Branch `feature/v1.3.0-rc5` at commit `afdcc821`, module `src/polyzymd/analyses/` (104 files, 51,982 lines, 1,222 tests). Reviewed 2026-09-11 by six independent adversarial review agents (statistics, observables, architecture, agent-facing API, docs/references/tests, loading/provenance) with cross-verification of the headline claims by the lead. Ground truth for statistical practice is Grossfield, Patrone, Roe, Schultz, Siderius and Zuckerman, "Best Practices for Quantifying Sampling Quality and Uncertainty in Molecular Simulations", LiveCoMS 1(1):5067 (2018), doi:10.33011/livecoms.1.1.5067. The six full reviews with file and line evidence are in `reviews/`.

Paths below are relative to `src/polyzymd/analyses/` unless they start with `src/`, `docs/` or `tests/`.

## 1. Verdict

Your memory is right about the plugins and wrong about the framework. The nine analysis packages (31,141 lines) are scripts fitted into a plugin contract. Each one carries its own aggregation, comparison, statistics, formatting and plotting, and every one bypasses the framework's generic paths. The framework itself (11,432 lines across `base.py`, `_framework/`, `mda/`, `orchestrator.py`, `stats.py`) is not script-like. It is over-built in the opposite direction, with two parallel lifecycles that define the same concepts twice.

The hard scientific decision is correct everywhere. Every cross-condition test and every condition-level error bar treats the replicate as the sampling unit. No per-frame hypothesis test exists. Equilibration is one global value, applied uniformly to every condition and replicate, and recorded in provenance. There is no path by which a diagnostic selects data. That satisfies the two rules in the LiveCoMS checklist that most often go wrong in the literature.

Three confirmed defects would not survive a referee who checks the numbers.

1. The correlation-time estimator every plugin calls reports a statistical inefficiency of 3.0 on white noise. I reproduced this: on 5,000 uncorrelated samples it reports 1,666 independent frames. The correct Chodera 2007 implementation sits in the same file and nothing calls it.
2. The hydrogen bond analysis hands the union of the protein and polymer selections to MDAnalysis as both the donor and the acceptor selection with no element restriction. Every carbon bonded to a hydrogen becomes a donor and every atom becomes an acceptor. I confirmed this in `hydrogen_bonds/_mda.py:159-166` and the plan builder at `:852-880`.
3. No confidence interval exists anywhere in the module. Every error bar is one standard error of the mean, unlabeled. With three replicates the 95 percent coverage factor is 4.30, so a reader who treats non-overlapping bands as significance is wrong by a factor of four.

The agent-facing problem is structural rather than a bug. There is no `polyzymd analyze` command, no `--format agent` for analyses, results carry no units except in RMSF, and reaching a first validated number from zero knowledge cost a simulated agent about 55,000 tokens of reading plus a 14 kB comparison YAML to edit. The `polyzymd status --format agent` command and its skill are the right template and were never applied here.

## 2. Findings verified by the lead

The reviewers each traced code. I independently reran the two most consequential claims and grepped three more.

| Claim | How verified | Result |
|---|---|---|
| `estimate_correlation_time(method="integration")` gives g = 3 on iid data | Ran in `pixi run -e analysis`, 30 series of 5,000 normal samples | g = 3.00, n_independent = 1,666 of 5,000. `statistical_inefficiency()` on the same data gives 1.03. On an AR(1) series with true g = 19 the plugin path gives 22.0 and the unused function 21.0. Bias is roughly +2 in g, dominated by the `max(tau, dt)` floor and the lag-0 half term. |
| Only the biased estimator is called | grep | Six callers (rmsd, rg, sasa, distances, catalytic_triad, rmsf). Zero callers of `statistical_inefficiency()`. |
| H-bond donors and acceptors are the whole group union | Read `hydrogen_bonds/_mda.py:159-166, 852-880`, MDAnalysis 2.10 `_get_dh_pairs` source | Confirmed. No `element N O` or `guess_donors` anywhere in the package. |
| Distances plugin applies no multiple-comparison correction | grep `fdr|adjust|benjamini` in `distances/__init__.py` | Zero hits. RMSD has nine. |
| Contacts aggregator claims autocorrelation correction it does not compute | Read `contacts/_aggregator.py:5-20`, grep | Docstring promises "autocorrelation-corrected uncertainties via statistical inefficiency" and cites Chodera 2007. The only statistics call in the file is a plain across-replicate SEM. |

A note on direction. The estimator bias overstates uncertainty for the within-replicate `sem_*` fields, so it is conservative. The comparison pipeline does not use those fields, so p-values are unaffected. The damage is that every `n_independent_frames` number and every "low statistical reliability" warning in the artifacts is wrong, and RMSF uses the inflated correlation time to discard at least half of its frames.

## 3. Compliance with the LiveCoMS checklist

| Checklist item | Status | Evidence |
|---|---|---|
| Do not cherry-pick; apply sampling metrics uniformly | Met | One global equilibration time (`src/polyzymd/config/comparison.py:58`), resolved once (`_framework/lifecycle.py:255, 578`), aggregation refuses replicates with mismatched windows (`mda/aggregation.py:410-439`). Convergence and N_eff never select data. |
| Remove an equilibration portion | Met, with a gap | Fixed user offset, default 10 ns. No automated detection (Chodera 2016) and no sensitivity report showing how results depend on the offset, which the article asks for. |
| Multiple runs; compare results across them | Met | Replicate is the unit in every test. Individual replicate points are drawn on bar plots for seven of nine plugins (not contacts or distances). |
| Estimate number of independent samples | Implemented incorrectly | Section 2, row 1. Block averaging is cited in `shared/autocorrelation.py:44-48` but does not exist. |
| Quantify uncertainty with confidence intervals | Not met | No `t.ppf`, no 1.96, no bootstrap, no jackknife in the module. Only SEM. |
| Describe what error bars mean in figures and tables | Not met | No `errorbar` or `fill_between` call carries a label; no JSON key states the uncertainty kind. The CLI table header does say "SEM". |
| Report the full procedure, publish scripts and data | Partly met | Artifacts record equilibration, test, correction, alpha, replicates, file identities. Missing: per-plugin algorithm version, package versions in most artifacts, git commit. |
| Semiquantitative checks: time series, two halves | Partly met | Time series plots exist. The sliding-window convergence heuristic in `shared/convergence.py` is uncited ("adapted from a collaborator notebook"), runs only for RMSD, and its default slope threshold of 5e-4 Å/ns sits below the noise of successive window means, so the flag tracks noise. |

Additional statistical findings, all confirmed:

- `ttest_method: welch` is silently ignored by rmsd, rg, sasa, contacts and distances. Only `stats.py:181` threads it through. With n = 3 and heteroscedastic conditions this changes p-values.
- The Benjamini-Hochberg family is defined differently across code paths (all metrics pooled in `stats.py:509-523`, per-analysis in contacts, absent in distances). The ANOVA is described as an omnibus gate but gates nothing.
- Cohen's d is reported with adjective labels at n = 3 without the Hedges small-sample correction.
- Singleton conditions write `sem = 0.0` to JSON rather than null.
- RMSF in external-reference mode computes per-atom RMSD about a crystal structure and stores it in the same `rmsf_values` field as true fluctuations (`rmsf/_mda.py:454-465`). The docs explain the distinction; the data format does not.

## 4. Observable-level correctness

Ordered by consequence for a polymer-enzyme paper. Full detail in `reviews/review_observables.md`.

1. Hydrogen bonds admit carbon donors and acceptors (high). PEG or polyacrylate with no polar hydrogens still produces "polymer donor" events. Because the fraction of such geometries differs by polymer chemistry, cross-condition comparisons are biased rather than offset. Fix by restricting both selections to `(<union>) and (element N O)` (add S if wanted), keeping `element H` for hydrogens, and recording the effective selections in provenance. Cite Smith et al. 2019 for the MDAnalysis implementation and Arunan et al. 2011 for the IUPAC definition.
2. Truncated or in-progress final segments load and cache silently (high). `src/polyzymd/engines/openmm/engine.py:263-270` includes any non-empty segment DCD regardless of status, and `_framework/lifecycle.py:376` reuses a cached replicate result with no size or mtime check. An analysis started mid-campaign freezes at a partial window and `finalize` and `plot-all` keep reusing it.
3. Engine-dependent periodic boundary semantics (high). The OpenMM path applies no unwrap or make-whole, while the GROMACS path prefers `prod_centered.xtc` produced with `-pbc mol -center`. The same plugin therefore sees different coordinates by engine, recorded only in a filename. Rg in fragment mode and long pair distances are exposed if a molecule straddles the boundary.
4. Bonds are silently absent for PDB topologies above 99,999 atoms (medium, suspected for typical solvated systems). OpenMM writes hex serials, MDAnalysis then refuses all CONECT records. Rg fragment mode falls back to one fragment and contacts assigns all polymer residues to chain 0, each with only a warning.
5. Distances and triad align the trajectory in memory, then apply minimum image with the unrotated box (medium, impact conditional). Harmless for short triad distances, wrong for long centre-of-mass pairs. Alignment is unnecessary for distances and copies the entire solvated system into memory.
6. RMSD default reference is a per-replicate centroid (medium, design choice). Each replicate is measured against its own structure, so cross-condition RMSD compares spread, not deviation from a common state. The docs recommend `external` for comparisons; the default and the metric name do not reveal the mode. The centroid search is also a single Kabsch pass to the start frame and ignores the analysis stride.
7. DSSP counts residues mdtraj cannot assign as coil (low-medium). Reachable with Amber or CHARMM residue names and modified conjugation-site residues.
8. Contacts apply the 4.5 Å heavy-atom convention to all atoms including hydrogens, undocumented (low-medium).
9. Empty selections yield zero results instead of errors in contacts, hydrogen bonds and SASA (low-medium). Rg does this correctly with an explicit skipped payload.
10. In-memory alignment shifts the time axis by one frame relative to plugins that read DCD timestamps (low).

Done well at this level: the segment lineage check refuses overlapping or gapped restart segments, equilibration resolves in absolute trajectory time with the DCD first-frame offset handled, units are consistent (Å throughout, explicit nm conversions for mdtraj SASA), contacts and pair distances use minimum image with the per-frame box, SASA uses documented mdtraj defaults with solvent excluded, the maximum-ASA table cites Tien et al. 2013, and each plugin gets a fresh universe so in-place alignment never leaks.

## 5. Architecture and extensibility

Full detail with counts in `reviews/review_architecture.md`.

Two frameworks define the same concepts. `Analysis` in `base.py` with `_framework/` and the `mda/` artifact layer each own a lifecycle, aggregation, comparison and contexts. The module has 8 `*Context` classes, 23 `*Result`, 15 `*Summary`, 25 `*Settings`, 17 `*Error` and 11 `*Collector` classes, three lifecycle modules, two `ANOVAResult` classes, and a `default_compare` that branches on which result family it received. More than 30 "stale cache" error strings exist to police the seam.

Statistics are done per plugin, not once. There are seven comparison engines. The hooks `aggregate`, `build_mda_jobs`, `build_mda_collector`, `compare`, `plot` and `format` are overridden by nine of nine plugins. The framework's generic `aggregate` is used by none. This is where the Welch, FDR and N_eff inconsistencies in section 3 come from: whether a comparison is corrected depends on which plugin you ran.

Contributor burden is about ten times MDAnalysis. Running `polyzymd new-analysis probe --dry-run` generates a 337-line plugin (three classes, nine functions) and a 277-line test. The contributor hand-builds provenance inside a collector. Helpers such as `mdanalysis_version` (seven copies), `_combined_warnings` (ten copies) and `_is_boolean_frame_value` (four copies) are pasted rather than imported. MDAnalysis needs about 30 lines for the same analysis.

Duplication across plugins runs 15 to 33 percent by file role (token-normalised 5-gram Jaccard). A 101-line `plot_rmsd_comparison_bars` is byte-identical to `plot_rg_comparison_bars` apart from the metric name. `_plot_settings.py` for rmsd and sasa have Jaccard 1.00. `_apply_fdr_correction` is copied verbatim into three `__init__.py` files although `shared/multi_run_comparison.apply_fdr_correction` exists.

Layering has cycles: `base` and `_framework`, `base` to `mda` to `_framework`, `mda` to `stats` to `base`, `shared` to `mda`, and `config` to `analyses` (parsing a comparison YAML imports every plugin). `base.py` rewrites `__module__` on 15 classes. `Analysis` inherits `ABC` with no abstract methods; the contract is enforced by a `TypeError` at class definition whose main job is rejecting hooks removed in earlier refactors. There are 326 `: Any` parameters and 115 `-> Any` returns.

Out-of-tree extension is impossible. Discovery scans `polyzymd.analyses.__path__` only, there is no entry-point group, and the scaffold writes into the repository. The documentation calls analyses "the primary extensibility axis".

Config coupling is narrow in practice but mandatory in types. `TrajectoryLoader` uses three `SimulationConfig` fields, but `Condition`, `ReplicateContext` and the loader all require the full object. Nothing accepts a topology path plus trajectory paths. That is why external agents only use the universe loader and then write their own loop.

Caching is inverted. Aggregates and comparisons are reused by `finalize`, `plot-all` and the HPC worker commands with no version or trajectory-identity check, although both are recorded. Meanwhile `compare run` always recomputes every replicate; no plugin reads `ctx.recompute`. The docs claim the opposite.

Done well: all eight `except Exception` sites wrap and re-raise typed errors, artifacts are pydantic models with a schema version and SHA-256 sidecars, `FrameSelection` and `MDABackendPolicy` are sound, discovery is deterministic, matplotlib is lazy in all plugins (seaborn is not, see section 6), and the framework core has a test-to-source ratio near 0.8.

## 6. Agent-facing API

Full detail and the token log in `reviews/review_api_agents.md`.

A reviewer played an agent from zero knowledge, using only what the repository exposes. Reaching four goals (load a universe, run RMSF with replicate uncertainty for one condition, compare Rg across two conditions as JSON, add a new analysis) took about 25 tool calls and 55,000 to 60,000 tokens of reading. The first three attempts failed because the default pixi environment has no `polyzymd` and nothing but a README table says to use `-e analysis`. The natural guess `polyzymd analyze` does not exist. A single-condition RMSF requires writing a full comparison YAML.

Two public universe loaders exist with no written preference. `TrajectoryLoader.load_universe` in `shared/loader.py` is the de facto canonical one; `UniverseProvider.load_universe` in `mda/universe.py` wraps it and adds provenance. The contributor guide recommends one and the architecture page presents the other.

The JSON output is a good envelope (schema version, package versions, config hash, source files, warnings, adjusted p, effect size, direction, a `testable` flag) but an agent cannot answer "by how much, in what units, with what interval, from how many frames" without reading code. The only `"unit"` key in the module is in RMSF. `MetricValue` has no unit field.

`list_analyses()` takes 2.9 seconds because `catalytic_triad/_plotters.py:34` imports seaborn at module level, four lines below a comment saying seaborn is lazy.

Proposed shape, reusing `run_comparison`, `ComparisonConfig` and `PairwiseResult`:

```
polyzymd analyze <name> -c A/config.yaml [-c B/config.yaml ...] [--eq 10ns] [--format agent|json|table]
```

One `-c` gives a per-condition summary; two or more give pairwise comparisons with the first as control. A Python equivalent `polyzymd.analyses.protocols.analyze(name, configs, ...)` returns a `ProtocolReport` whose every number carries `unit`, `ci95`, `ci_method`, `test`, `correction`, `frames_per_replicate`, provenance and a one-line `verdict`. The `--format agent` renderer prints at most 25 lines. A 40-line `SKILL.md` mirroring `polyzymd-status` tells the agent the one command, the environment, the verdict vocabulary, and "do not write your own MDAnalysis loop". Estimated cost per analysis drops from roughly 55,000 tokens to 1,000 to 2,000. The only new statistics are a t-based interval and a unit field.

## 7. Loading and provenance

Full detail in `reviews/review_loading_provenance.md`.

What "load from config" resolves to: `solvated_system.pdb` written at build time from the OpenFF topology (pre-minimisation coordinates, elements present, bonds only via CONECT for non-standard residues), plus every `production_N/production_N_trajectory.dcd` ordered by integer index with contiguity enforced. Time comes from the DCD header only. The loader wraps the ChainReader in `_TimestampPreservingTrajectory` to keep absolute times; this object is not an MDAnalysis `ProtoReader`, so type-checking library code may misbehave.

Provenance present in MDA-layer artifacts: file paths, sizes and mtimes for every input, full frame window with equilibration and time reference, settings fingerprint, selections in most plugins. Missing: per-plugin algorithm version, MDAnalysis, numpy and mdtraj versions in replicate artifacts (the manifest field exists but `write_manifest` is never called), git commit, input content hashes, alignment reference identity. Random seeds are not needed; the module has no RNG.

A result can be traced to exact frames of exact files, with indices referring to the ChainReader concatenation. It cannot be traced to the segment and local frame without reopening files.

## 8. Documentation, tests and attribution

Full detail with page-by-page evidence in `reviews/review_docs_refs_tests.md`.

Docs are accurate on imports, CLI flags and YAML keys (a delegated pass verified all 29 `api/analyses_mda` symbols and found zero dead links). Reference pages are stale in specifics: `--eq-time` documented as default `0ns` on five pages when the effective default is `10ns`; sidecar file names, payload keys and class names that do not exist (`DistancePair`, `TriadPair`, `timeseries_sidecar`, `profile_sidecar`, `hbonds_eq*.json`); triad `simultaneous_contact_fraction` shown as a fraction but stored as a percent; four quickstarts run `polyzymd info -c config.yaml --scratch-dir` when `info` takes no options.

Only one page defines what SEM means. No plugin page distinguishes across-replicate SEM from within-replicate `std/sqrt(n_independent)`, although the code produces both and plots mix them. Two explanation pages mention confidence intervals and block averages as if produced.

The contributor guide contradicts the scaffold it documents (`METRIC_NAME = "mean_shell_count"` versus the generated `solvent_shell_value`), omits `extract_metrics()`, and never asks a contributor to cite the method, write a known-answer test, or state what the error bar means.

Tests: 1,222 collected, over 95 percent plumbing. Real known-answer tests exist for RMSD of a stretched toy, Benjamini-Hochberg, ANOVA regression, convergence plateau and window rounding. Zero tests for `statistical_inefficiency*`, `compute_sem` correctness, Kabsch, SASA geometry or DSSP (mocked). No real trajectory fixture; the GROMACS "smoke" tests pass `universe=object()`. No physical-validation test (two halves agree, SEM calibration). The first two missing tests in the reviewer's list would fail today against the estimator bug.

Attribution: no `CITATION.cff` and no references page. MDAnalysis, mdtraj, DSSP, Shrake-Rupley, Kabsch and QCP, Welch, Tukey, Cohen and SciPy are used uncited. `statistical_inefficiency()` is a structural port of pymbar's algorithm without acknowledgement (the `_multiple` variant does acknowledge it, but carries no MIT notice). `scikit-learn` is a declared dependency never imported by the module.

## 9. What to change, in order

The correctness fixes are days of work. The structural collapse is about seven weeks. My recommendation is that v1.3.0 ships correctness, the agent protocol, provenance and references, and that the framework collapse becomes v1.4.0 behind an unchanged CLI. Shipping wrong N_eff values and carbon hydrogen bonds for another release cycle costs more than shipping a framework that is larger than it should be.

### Phase 0, correctness (2 to 3 days, all in v1.3.0)

1. Make `estimate_correlation_time` a thin wrapper around `statistical_inefficiency()`, remove the `max(tau, dt)` floor, and re-pin tests to g ≈ 1 on white noise and g ≈ (1+φ)/(1−φ) on AR(1).
2. Stop subsampling RMSF frames; keep τ as a diagnostic; put uncertainty on RMSF via across-replicate SEM.
3. Restrict hydrogen bond donors and acceptors to N and O (S optional) and record the effective selections.
4. Thread `ttest_method` and `posthoc_method` through every plugin `compare`; add FDR to distances; define the BH family once and apply it to ANOVA consistently.
5. Add `ci95_low`, `ci95_high`, `ci_method` and `unit` to every metric model; label every error bar and band ("mean ± 1 SEM over n replicates, production window t ≥ X ns"); add an `uncertainty` block to every artifact.
6. Consult `progress.json` status per segment before loading; refuse cached replicates whose recorded size or mtime differ from disk.
7. Declare a `pbc` policy on load, record it, and fail (not warn) when fragment mode or chain identity needs bonds that are absent.
8. Drop alignment in distances and triad, or disable minimum image when aligned.
9. Fix the contacts aggregator docstring or implement what it claims; delete the Flyvbjerg citation until block averaging exists.

### Phase 1, agent protocol (3 to 4 days, v1.3.0)

`polyzymd analyze`, `protocols.analyze()`, `ProtocolReport`, `--format agent`, the skill file, and a `load(config, replicate, window=...)` plus `load_files(topology, trajectories, dt_ps=...)` pair returning a `LoadedTrajectory` with a `LoadProvenance` object. Write down in one place which loader is canonical. Fix the seaborn import.

### Phase 2, trust (1 week, v1.3.0)

Known-answer tests: white noise and AR(1) for both estimators; SEM calibration over synthetic replicates; Kabsch on a rotated cloud; single-sphere SASA against 4π(r+p)²; DSSP on an ideal helix; Rg on a cube; contacts on a two-atom toy; one committed tiny real trajectory (≤ 50 frames, ≤ 500 atoms) exercising equilibration and stride end to end and cross-checking RMSF against `MDAnalysis.analysis.rms.RMSF`. Add `CITATION.cff`, `docs/source/explanation/references.md`, and a NumPy-style `References` section plus a one-line "Method" statement in every plugin `__init__.py`, printed under the results table by the default formatter. Fix the stale reference pages listed in section 8.

### Phases 3 to 5, framework collapse (about 6 weeks, v1.4.0)

Introduce `Observable(name, unit, kind, values, index)` where `kind` is one of `mean_of_timeseries`, `fluctuation`, `fraction`, `distribution`, `profile`, and reduce a plugin to `Settings` plus `compute(universe, frames, settings) -> Sequence[Observable]`. The framework then owns, once: persistence with a framework-written identity block (polyzymd version, plugin source hash, settings fingerprint, config hash, input file identity), aggregation by kind (N_eff-corrected mean and SEM, fluctuation, binomial fraction, KDE, per-index profile), comparison (Welch or Student plus BH or Tukey from one defaults object), replicate-cache reuse keyed on the identity block, and generic plotting and formatting per kind with an optional `extra_plots` hook. Collapse the `MDA*Context` classes into the four framework contexts, delete the `ConditionSummary`/`ComparisonResult` family in favour of the artifact envelope, remove the `__module__` rewriting and the removed-hook police, give `Analysis` real abstract methods or make it a Protocol, add an entry-point group `polyzymd.analyses` for out-of-tree plugins, and replace `SimulationConfig` in contexts with a `TrajectorySource` protocol. Port rmsf, rg, rmsd and sasa first (all scalar per run today), then plotters, then contacts, hydrogen bonds and triad with an optional `Aggregator` protocol for residence times.

Expected outcome: roughly 15,000 to 17,000 of the 52,000 lines disappear, the new-plugin floor drops from about 340 scaffold lines to about 40, and statistical practice becomes a framework guarantee rather than a per-plugin habit. Keep untouched: `mda/artifacts.py`, `mda/store.py`, `mda/frame_selection.py`, `mda/job.py`, `mda/universe.py`, `shared/inferential_statistics.py`, `shared/autocorrelation.py` (after the fix), `shared/loader.py` internals, `discovery.py`, `exceptions.py`.

## 10. References to add

Methods the module implements or adapts, with what it currently cites. Entries marked "keep" are already cited correctly.

| Method | Where used | Cite |
|---|---|---|
| MDAnalysis | loader, alignment, RMSD, RMSF, contacts, H-bonds | Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787; Gowers et al. 2016, Proc SciPy, doi:10.25080/Majora-629e541a-00e |
| MDAnalysis HydrogenBondAnalysis | `hydrogen_bonds/_mda.py` | Smith et al. 2019, PCCP 21:9845, doi:10.1039/C9CP01532A |
| IUPAC hydrogen bond definition | H-bond criteria | Arunan et al. 2011, Pure Appl Chem 83:1637, doi:10.1351/PAC-REC-10-01-02 |
| mdtraj (SASA, DSSP) | `sasa/_mda.py`, `secondary_structure/_mda.py`, `rmsf/_mda.py` | McGibbon et al. 2015, Biophys J 109:1528, doi:10.1016/j.bpj.2015.08.015 |
| Shrake-Rupley SASA and radii | `sasa/_mda.py` | Shrake and Rupley 1973, J Mol Biol 79:351, doi:10.1016/0022-2836(73)90011-9; Bondi 1964, J Phys Chem 68:441 |
| Maximum ASA | `shared/aa_classification.py` | Tien et al. 2013 (keep) |
| DSSP | `secondary_structure/_mda.py` | Kabsch and Sander 1983, Biopolymers 22:2577, doi:10.1002/bip.360221211 |
| Kabsch superposition | `shared/centroid.py` | Kabsch 1976, Acta Cryst A32:922, doi:10.1107/S0567739476001873 |
| QCP alignment (via MDAnalysis) | `shared/alignment.py`, `rmsd/_mda.py` | Theobald 2005, Acta Cryst A61:478, doi:10.1107/S0108767305015266; Liu, Agrafiotis and Theobald 2010, J Comput Chem 31:1561 |
| Statistical inefficiency, N_eff | `shared/autocorrelation.py` | Chodera et al. 2007, JCTC 3:26, doi:10.1021/ct0502864 (keep); Shirts and Chodera 2008, J Chem Phys 129:124105, doi:10.1063/1.2978177 plus pymbar MIT notice; Janke 2002, NIC Series 10:423; Sokal 1997, doi:10.1007/978-1-4899-0319-8_6 |
| Best practices, coverage factors | everywhere uncertainty is reported | Grossfield et al. 2018, LiveCoMS 1:5067 (keep, extend to code); Grossfield and Zuckerman 2009, Annu Rep Comput Chem 5:23, doi:10.1016/S1574-1400(09)00502-7; JCGM 100:2008 (GUM) |
| Equilibration detection (if added) | convergence, window | Chodera 2016, JCTC 12:1799, doi:10.1021/acs.jctc.5b00784; Yang et al. 2004; Klimovich et al. 2015 |
| Block averaging (only if implemented) | none today | Flyvbjerg and Petersen 1989, J Chem Phys 91:461, doi:10.1063/1.457480 |
| Bootstrap (only if implemented) | none today | Efron 1979, Ann Stat 7:1, doi:10.1214/aos/1176344552; Efron and Tibshirani 1993 |
| Structural decorrelation, effective sample size | RMSF subsampling rationale | Lyman and Zuckerman 2007, J Phys Chem B 111:12876, doi:10.1021/jp073061t; Zhang, Bhatt and Zuckerman 2010, JCTC 6:3048, doi:10.1021/ct1002384 |
| Welch t-test | `shared/inferential_statistics.py` | Welch 1947, Biometrika 34:28, doi:10.1093/biomet/34.1-2.28 |
| Tukey HSD | `shared/inferential_statistics.py` | Tukey 1949, Biometrics 5:99, doi:10.2307/3001913 |
| Benjamini-Hochberg | `shared/inferential_statistics.py` | Benjamini and Hochberg 1995 (keep) |
| Effect size | `shared/inferential_statistics.py` | Cohen 1988, Statistical Power Analysis for the Behavioral Sciences, 2nd ed.; Hedges 1981, J Educ Stat 6:107, doi:10.3102/10769986006002107 |
| SciPy | all tests | Virtanen et al. 2020, Nat Methods 17:261, doi:10.1038/s41592-019-0686-2 |
| Gaussian KDE | `distances/_mda.py` | Scott 1992, Multivariate Density Estimation, Wiley |
| Convergence heuristic | `shared/convergence.py` | Name the collaborator or replace with a cited method (Hess 2002, Phys Rev E 65:031910; Knapp et al. 2018 already cited in RMSF docs) |
| RMSF vs B-factor | `rmsf/` docs | Kuzmanic and Zagrovic 2010, Biophys J 98:861 |

## 11. Method and limits

Each reviewer read code rather than docstrings and marked findings as confirmed (traced or executed) or suspected. I reran the estimator claim and the hydrogen bond claim myself and grepped three others. No real trajectory exists on this machine, so nothing was executed against production data; findings about large-system bond parsing, restart-boundary duplicates and periodic boundary splitting are inferred from the OpenMM and MDAnalysis sources and marked suspected where the reviewer could not exercise them. Duplication percentages are token-normalised estimates. One reviewer was killed by a rate limit and relaunched; its report is complete. The repository was not modified.

## 12. Executing the plan with Claude Code

Sources: the Claude Code best-practices page (code.claude.com/docs/en/best-practices), the `/goal` page (code.claude.com/docs/en/goal), the hooks guide, the workflows page, and the routines page. Local install is Claude Code 2.1.269, which has `/goal`.

The unit of work is one session, one branch, one stacked PR layer. `/goal` is the right tool inside a session and the wrong tool across sessions. Its condition is judged by a small evaluator model after each turn, it survives only `--resume`, it clears on context overflow or auth errors, and it cannot open a second PR. So the multi-week plan lives in a checklist file in the repo, and each session picks one item, sets a `/goal` whose condition is machine-checkable, and ends with a clean commit on its own branch.

Do these four things before the first fix session.

1. Make the repo legible to a fresh session. `.opencode/instructions/analysis-module.md` lists files that no longer exist (`_results.py`, `_cache.py`, `_paths.py`, `_plotting.py`). Rewrite it from the real tree, add a `CLAUDE.md` that imports `AGENTS.md` with an `@AGENTS.md` line so Claude Code reads the same rules Codex does, and state the two pixi environments (`analysis` for running, `test` for pytest) since the default environment has no numpy and cost the API reviewer three failed commands.
2. Write the known-answer tests first, as failing tests. This is the best-practices "write a failing test that reproduces the issue" pattern and it turns every Phase 0 item into a `/goal` condition. Put them in `tests/analyses/scientific/`: white noise and AR(1) for both correlation estimators, SEM calibration, Kabsch on a rotated cloud, single-sphere SASA, DSSP on an ideal helix, Rg on a cube, contacts on a two-atom toy, hydrogen bonds on a toy with one N-H···O and one C-H···O pair (exactly one bond expected), and Welch versus Student on unequal-variance fixtures for each plugin. Two of these fail today against the estimator bug, one against the hydrogen bond bug.
3. Add a deterministic gate. A `PreToolUse` hook on `Bash` that blocks `git commit` unless `ruff check src`, `black --check src` and the fast scientific test subset exit 0. Hooks are the documented way to make a check non-negotiable; a `/goal` condition is advisory and a Stop hook can be overridden after eight blocks.
4. Encode the standard as a repo skill. `.claude/skills/livecoms-check/SKILL.md` holding the LiveCoMS checklist items from section 3, the "replicate is the unit" rule, the required `unit`, `ci95`, `ci_method` and `uncertainty` fields, and the citation convention from section 10. Every session invokes it before committing, and the reviewer subagent uses it as its rubric. A second skill, `polyzymd-analyze`, is the agent-facing protocol from section 6 and is itself a Phase 1 deliverable.

A session then looks like this. Start with `claude --worktree fix-correlation-time` on top of the stack's parent branch. Open in plan mode, read the checklist item and the audit section, and write the plan. Set the goal:

```
/goal `pixi run -e test pytest tests/analyses/scientific -q -k "correlation or inefficiency"` exits 0, `ruff check src` and `black --check src` are clean, all pre-existing tests in tests/analyses/shared still pass, and the change is committed on this branch with a conventional-commit subject under 50 characters and no attribution trailer
```

Then let it run. When it stops, run `/code-review` so a fresh-context reviewer sees only the diff and the rubric, fix what survives, check the item off in the checklist file, and stop. Adopt the branch into the stack with `gh stack`, and do not open the PR until you have confirmed completion, per the repo rules.

Parallelism is worth it only for independent mechanical items. Threading `ttest_method` through five plugins, adding `References` sections to nine `__init__.py` files, and fixing the stale reference pages in section 8 are each many-files-same-change jobs. Those suit a workflow or `/batch`: one subagent per plugin, each in its own worktree, each producing a small PR, then a single verify pass. The framework collapse in v1.4 is the opposite. It is design-heavy, touches shared files, and each plugin port depends on the `Observable` contract landing first. Run it sequentially: one plan-mode session to write the design document and the `Observable` contract with its tests, then one session per plugin port, rmsf first because it is smallest and already declares a statistical policy.

Two optional additions. A cloud routine (`/schedule`, research preview) can run the full slow test suite and the docs build nightly against the stack tip and open an issue on failure. Agent teams and agent view are experimental and off by default; nothing here needs them.

Session count, rough: Phase 0 in six sessions (estimator plus RMSF subsampling; hydrogen bonds; Welch and FDR threading via fan-out; CI fields and labels; segment status and cache freshness; PBC policy and alignment removal). Phase 1 in three (loader pair; `analyze` command and `ProtocolReport`; skill and docs). Phase 2 in four (tests are largely written by then; references page and `CITATION.cff`; stale docs via fan-out; contributor guide). v1.4 in fifteen to twenty. Each session should end with the repo mergeable and the checklist updated, so a rate-limit death like the one that killed a reviewer in this audit costs one item, not the plan.
