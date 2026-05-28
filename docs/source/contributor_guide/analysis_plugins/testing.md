# Test your analysis plugin contribution

This how-to shows the focused tests to write before opening an analysis plugin
pull request. Use it after you have a working plugin shape and need confidence
that discovery, settings, MDAnalysis jobs, artifacts, aggregation, comparison,
sidecars, and plots follow the PolyzyMD contract.

The snippets use a teaching plugin named `solvent_shell`. Replace that name with
your plugin name and keep the tests in `tests/analyses/plugins/test_<name>.py`.

## Run the focused test file

Run one plugin test file while you iterate:

```bash
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
```

Before a pull request, also run the broader analysis test subset and formatting
checks listed in {doc}`checklist`.

## Test discovery and settings validation

Discovery tests catch missing class variables, naming mistakes, and plugin files
that cannot be imported. Settings tests should cover defaults and at least one
invalid value.

```python
import pytest

from polyzymd.analyses.discovery import clear_cache, get_analysis, list_analyses
from polyzymd.analyses.solvent_shell import SolventShellAnalysis, SolventShellSettings


def test_solvent_shell_is_discoverable() -> None:
    clear_cache()

    assert list_analyses()["solvent_shell"] is SolventShellAnalysis
    assert get_analysis("solvent_shell") is SolventShellAnalysis


def test_settings_defaults_and_validation() -> None:
    settings = SolventShellSettings()
    assert settings.selection == "protein and name CA"
    assert settings.scale == 1.0

    with pytest.raises(ValueError):
        SolventShellSettings(selection="")

    with pytest.raises(ValueError):
        SolventShellSettings(scale=0.0)
```

Keep discovery tests lightweight. They should not load trajectories or require
external simulation data.

## Test `build_mda_jobs()` with fakes

Use small fake universe objects or mocks so the test verifies job construction,
frame selection, and settings wiring without importing real trajectories.

```python
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from polyzymd.analyses.base import Condition
from polyzymd.analyses.mda import FrameSelection, MDAUniversePolicy
from polyzymd.analyses.solvent_shell import SolventShellAnalysis, SolventShellSettings


class FakeTrajectory:
    def __len__(self) -> int:
        return 10


class FakeUniverse:
    trajectory = FakeTrajectory()

    def select_atoms(self, selection: str):
        return SimpleNamespace(n_atoms=3, positions=[])


def test_build_mda_jobs_returns_named_jobs() -> None:
    condition = Condition(
        label="Control",
        config_path=Path("/fake/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    ctx = SimpleNamespace(
        universe=FakeUniverse(),
        settings=SolventShellSettings(selection="protein and name CA", scale=1.0),
        frame_selection=FrameSelection(start=0, stop=10, step=2, n_frames_total=10),
        replicate_context=SimpleNamespace(condition=condition),
        replicate=1,
        universe_policy=MDAUniversePolicy(condition_label="Control", replicate=1),
    )

    jobs = SolventShellAnalysis().build_mda_jobs(ctx)

    assert jobs
    assert jobs[0].name == "solvent_shell"
    assert jobs[0].frame_selection.step == 2
```

If your job uses an `AnalysisBase`-compatible worker, patch that worker with a
small fake class. If your job uses `MDAAnalysisJob.from_function()`, assert that
the function kwargs reflect `ctx.settings` and that explicit frame selections are
preserved.

## Test that the collector returns a `ReplicateArtifact`

Collectors translate runtime job output into PolyzyMD's durable artifact
contract. The test should assert the artifact type, scalar metrics, provenance,
warnings, and any registered sidecars.

```python
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from polyzymd.analyses.base import Condition
from polyzymd.analyses.mda import (
    ArtifactStore,
    FrameSelection,
    MDABackendPolicy,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.solvent_shell import SolventShellAnalysis, SolventShellSettings


def test_collector_returns_replicate_artifact(tmp_path: Path) -> None:
    condition = Condition(
        label="Control",
        config_path=Path("/fake/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    replicate_context = SimpleNamespace(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        result_path=tmp_path / "run_1" / "result.json",
        settings=SolventShellSettings(),
    )
    replicate_context.output_dir.mkdir(parents=True)

    collector_ctx = MDACollectorContext(
        analysis_name="solvent_shell",
        replicate_context=replicate_context,
        frame_selection=FrameSelection(start=0, stop=5, step=1, n_frames_total=5),
        universe_policy=MDAUniversePolicy(condition_label="Control", replicate=1),
        artifact_store=ArtifactStore(replicate_context.output_dir),
        settings_fingerprint="test-settings",
    )
    job = MDAJobResult(
        name="solvent_shell",
        analysis=SimpleNamespace(),
        results={
            "metrics": {"mean_shell_count": 4.0},
            "mean_shell_count": 4.0,
            "n_frames": 5,
        },
        run_kwargs={},
        frame_selection=collector_ctx.frame_selection,
        backend_policy=MDABackendPolicy(),
        universe_policy=collector_ctx.universe_policy,
    )

    collector = SolventShellAnalysis().build_mda_collector(collector_ctx)
    artifact = collector(collector_ctx, [job])

    assert isinstance(artifact, ReplicateArtifact)
    assert artifact.analysis_name == "solvent_shell"
    assert artifact.payload["metrics"]["mean_shell_count"] == 4.0
    assert artifact.provenance
```

Do not assert that a collector serializes raw MDAnalysis `Results` objects. A
collector should map runtime containers into JSON-compatible payload values and
registered sidecars.

## Test aggregation and default metrics

For simple scalar plugins, write replicate artifacts to an `ArtifactStore`, then
call the plugin's aggregation path or default aggregation helper exposed by the
plugin. The key assertion is that the condition-level output contains a metric
summary with replicate values.

```python
from pathlib import Path
from unittest.mock import MagicMock

from polyzymd.analyses.base import AggregateContext, Condition, MetricValue
from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.solvent_shell import SolventShellAnalysis, SolventShellSettings


def test_aggregate_or_default_metrics(tmp_path: Path) -> None:
    analysis = SolventShellAnalysis()
    settings = SolventShellSettings()
    condition = Condition(
        label="Control",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )
    ctx = AggregateContext(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        equilibration="0ns",
        settings=settings,
    )
    frame_selection = {"start": 0, "stop": None, "step": 1, "frames": None}
    settings_fingerprint = analysis.aggregate_settings_fingerprint(settings)

    # The default Analysis.aggregate() path loads replicate artifacts from disk
    # under ctx.output_dir.parent / f"run_{replicate}" / "result.json".
    artifacts = [
        ReplicateArtifact(
            analysis_name="solvent_shell",
            condition_label="Control",
            replicate=1,
            payload={"metrics": {"mean_shell_count": 4.0}},
            provenance={"frame_selection": frame_selection},
            metadata={"settings_fingerprint": settings_fingerprint},
        ),
        ReplicateArtifact(
            analysis_name="solvent_shell",
            condition_label="Control",
            replicate=2,
            payload={"metrics": {"mean_shell_count": 6.0}},
            provenance={"frame_selection": frame_selection},
            metadata={"settings_fingerprint": settings_fingerprint},
        ),
    ]
    for artifact in artifacts:
        run_dir = ctx.output_dir.parent / f"run_{artifact.replicate}"
        ArtifactStore(run_dir).write_replicate_result(artifact)

    result = analysis.aggregate(ctx, artifacts)

    assert isinstance(result, ConditionArtifact)
    metric = result.payload["metrics"]["mean_shell_count"]
    assert metric["values"] == [4.0, 6.0]
    assert metric["mean"] == 5.0


def test_extract_metrics_for_custom_aggregate_model() -> None:
    metrics = SolventShellAnalysis().extract_metrics(
        {"mean_shell_count": 5.0, "sem_shell_count": 0.5, "replicate_values": [4.0, 6.0]}
    )

    if metrics:
        assert isinstance(metrics["mean_shell_count"], MetricValue)
```

The second test is relevant only when your plugin returns a custom aggregate
model or dictionary and therefore implements `extract_metrics()`. For the
canonical `ReplicateArtifact` to `ConditionArtifact` path, default comparison can
read condition artifact metrics directly.

## Test artifact-only plotting

Plot tests should create cached artifacts or sidecars, call `plot()`, and assert
that figures are written. They should not create a universe, read trajectories,
or run MDAnalysis jobs.

The scaffolded `SolventShellAnalysis` returns no plots until you implement a
`plot()` hook. Use this test pattern only after adding an artifact-only hook that
loads condition artifacts or sidecars and returns at least one output path. A
minimal hook can live in the plugin itself; the test below defines it inline only
to make the expected behavior explicit.

```python
from pathlib import Path
from unittest.mock import MagicMock

from polyzymd.analyses.base import Condition, PlotContext
from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact
from polyzymd.analyses.solvent_shell import SolventShellAnalysis, SolventShellSettings


class PlottingSolventShellAnalysis(SolventShellAnalysis):
    def plot(self, ctx: PlotContext) -> list[Path]:
        import matplotlib.pyplot as plt

        labels: list[str] = []
        means: list[float] = []
        for label, analysis_dir in ctx.analysis_dirs.items():
            artifact = ArtifactStore(analysis_dir).read_condition_result("result.json")
            metric = artifact.payload["metrics"]["mean_shell_count"]
            labels.append(label)
            means.append(float(metric["mean"]))

        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        path = ctx.output_dir / "mean_shell_count.png"
        fig, ax = plt.subplots()
        ax.bar(labels, means)
        ax.set_ylabel("Mean shell count")
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        return [path]


def test_plot_reads_condition_artifacts_only(tmp_path: Path) -> None:
    condition = Condition(
        label="Control",
        config_path=Path("/fake/config.yaml"),
        replicates=(1,),
        sim_config=MagicMock(),
    )
    analysis_dir = tmp_path / "Control" / "solvent_shell" / "aggregated"
    analysis_dir.mkdir(parents=True)
    ArtifactStore(analysis_dir).write_condition_result(
        ConditionArtifact(
            analysis_name="solvent_shell",
            condition_label="Control",
            payload={
                "metrics": {
                    "mean_shell_count": {"mean": 5.0, "sem": 0.2, "values": [4.8, 5.2]}
                }
            },
        ),
        "result.json",
    )
    ctx = PlotContext(
        conditions=(condition,),
        analysis_dirs={"Control": analysis_dir},
        results_dir=tmp_path / "comparison" / "solvent_shell",
        output_dir=tmp_path / "plots",
        settings=SolventShellSettings(),
    )

    paths = PlottingSolventShellAnalysis().plot(ctx)

    assert paths
    assert all(path.exists() for path in paths)
```

Prefer asserting artifact-only inputs and outputs over patching
`MDAnalysis.Universe`, because that patch can be fragile when plugins import
MDAnalysis lazily or never import it during plotting. Add a targeted trajectory
loader patch only if your plugin has a known plotting regression to guard.

## Validate sidecars where relevant

Plugins that write NPZ, CSV, or other sidecars should test both registration and
loading through `ArtifactStore`. This catches stale files, missing metadata, and
absolute-path leaks.

```python
from pathlib import Path

import numpy as np

from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact


def test_npz_sidecar_is_registered_and_validated(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path)
    sidecar = store.write_npz_sidecar(
        "sidecars/shell_counts.npz",
        counts=np.asarray([4.0, 6.0], dtype=np.float64),
        metadata={"kind": "shell_counts", "units": "count"},
    )
    artifact = ReplicateArtifact(
        analysis_name="solvent_shell",
        condition_label="Control",
        replicate=1,
        payload={"n_frames": 2, "metrics": {"mean_shell_count": 5.0}},
        sidecars=[sidecar],
    )
    store.write_replicate_result(artifact, "result.json")

    loaded = store.read_replicate_result("result.json")
    assert loaded.sidecars[0].path == "sidecars/shell_counts.npz"
    assert loaded.sidecars[0].metadata["kind"] == "shell_counts"

    with store.load_npz_sidecar(loaded.sidecars[0]) as data:
        assert np.allclose(data["counts"], [4.0, 6.0])
```

Add a negative test when your sidecar contract is safety-critical, for example by
removing a required array key or checking that your loader rejects incompatible
metadata.

## Label slow tests clearly

Most plugin tests should use fakes, `ArtifactStore`, and small arrays. If a test
needs real trajectories or external analysis data, mark it as slow and document
the data requirement in the test name or docstring.

```python
import pytest


@pytest.mark.slow
def test_solvent_shell_real_trajectory_smoke() -> None:
    """Smoke test requiring external trajectory data."""
    ...
```

Keep slow tests separate from the contributor's fast unit-test loop.

## Finish with the PR checks

Run these before asking for review:

```bash
pixi run -e build pytest tests/analyses/plugins/test_solvent_shell.py -v
pixi run -e build pytest tests/analyses/ -v
pixi run -e build ruff check src/ tests/analyses/plugins/test_solvent_shell.py
pixi run -e build black src/ tests/analyses/plugins/test_solvent_shell.py --check
pixi run -e build make -C docs clean html
```

Use {doc}`checklist` as the final review pass.
