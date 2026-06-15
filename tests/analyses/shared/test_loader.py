"""Regression tests for TrajectoryLoader engine-aware refactor.

Tests cover:
- OpenMM daisy-chain and non-canonical directory layouts
- GROMACS flat production layout
- GRO topology chain-ID warning
- engine_override parameter
- Lazy engine creation
- find_topology() backward compatibility for direct callers
- _infer_replicate() edge cases
- Error paths (missing topology, missing trajectories, missing working dir)
- Engine resolution error propagation
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from polyzymd.analyses.shared.loader import TrajectoryLoader, enrich_universe_elements

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_openmm_config(tmp_path: Path, *, engine: str = "openmm") -> MagicMock:
    """Build a minimal mock SimulationConfig for OpenMM tests."""
    config = MagicMock()
    config.engine = engine
    config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
    config.output.effective_scratch_directory = tmp_path
    config.generate_system_name.return_value = "test_system"
    return config


def _make_gromacs_config(tmp_path: Path) -> MagicMock:
    """Build a minimal mock SimulationConfig for GROMACS tests."""
    config = MagicMock()
    config.engine = "gromacs"
    config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
    config.output.effective_scratch_directory = tmp_path
    config.generate_system_name.return_value = "test_system"
    return config


def _create_openmm_daisy_chain(run_dir: Path, n_segments: int = 3) -> None:
    """Populate *run_dir* with a daisy-chain OpenMM directory layout."""
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "solvated_system.pdb").write_text("ATOM solvated")
    for i in range(n_segments):
        seg_dir = run_dir / f"production_{i}"
        seg_dir.mkdir()
        (seg_dir / f"production_{i}_trajectory.dcd").write_bytes(b"DCD")
        (seg_dir / f"production_{i}_topology.pdb").write_text("ATOM topology")


def _create_openmm_noncanonical(run_dir: Path) -> None:
    """Populate *run_dir* with a non-canonical OpenMM single-directory layout."""
    run_dir.mkdir(parents=True, exist_ok=True)
    prod_dir = run_dir / "production"
    prod_dir.mkdir()
    (prod_dir / "production_topology.pdb").write_text("ATOM topology")
    (prod_dir / "production_trajectory.dcd").write_bytes(b"DCD")


def _create_gromacs_flat(run_dir: Path, *, use_gro: bool = False) -> None:
    """Populate *run_dir* with a flat GROMACS production layout."""
    gromacs_dir = run_dir / "gromacs"
    gromacs_dir.mkdir(parents=True, exist_ok=True)
    if use_gro:
        (gromacs_dir / "solvated_system.gro").write_text("GRO content")
    else:
        (gromacs_dir / "solvated_system.pdb").write_text("ATOM topology")
    (gromacs_dir / "prod.xtc").write_bytes(b"\x00")


class _FakeTrajectory:
    """Minimal trajectory double with frame and time metadata."""

    def __init__(
        self,
        times_ps: list[float],
        *,
        current_frame: int = 0,
        expose_time: bool = True,
        raw_times_ps: list[float] | None = None,
    ) -> None:
        self.times_ps = times_ps
        self.raw_times_ps = raw_times_ps
        self.frame = current_frame
        self.expose_time = expose_time
        self.ts = SimpleNamespace(frame=current_frame, data={})
        self._update_timestep()

    def __len__(self) -> int:
        """Return the number of fake frames."""

        return len(self.times_ps)

    def __getitem__(self, frame_index: int) -> SimpleNamespace:
        """Move to a fake frame and return a timestamp object."""

        if frame_index < 0 or frame_index >= len(self.times_ps):
            raise IndexError(frame_index)
        self.frame = frame_index
        self._update_timestep()
        return SimpleNamespace(time=self.times_ps[frame_index], frame=frame_index)

    @property
    def time(self) -> float:
        """Return the current fake frame time."""

        if not self.expose_time:
            raise AttributeError("time")
        return self.times_ps[self.frame]

    def _update_timestep(self) -> None:
        """Update raw timestep metadata for the current fake frame."""

        self.ts.frame = self.frame
        if self.raw_times_ps is None:
            self.ts.data = {}
            return
        self.ts.data = {"time": self.raw_times_ps[self.frame]}


class _ElementFakeAtoms:
    """Minimal atom container for element-enrichment tests."""

    def __init__(
        self,
        *,
        names: list[str],
        types: list[str] | None = None,
        resnames: list[str] | None = None,
        elements: list[str] | None = None,
    ) -> None:
        self.names = names
        self.types = types if types is not None else []
        self.resnames = resnames if resnames is not None else ["ALA"] * len(names)
        if elements is not None:
            self.elements = elements

    def __len__(self) -> int:
        """Return the number of fake atoms."""

        return len(self.names)


class _ElementFakeUniverse:
    """Minimal universe that records added element topology attributes."""

    def __init__(self, atoms: _ElementFakeAtoms) -> None:
        self.atoms = atoms
        self.added_topology_attrs: dict[str, list[str]] = {}

    def add_TopologyAttr(self, name: str, values: list[str]) -> None:
        """Record an added topology attribute on the fake atoms object."""

        self.added_topology_attrs[name] = values
        setattr(self.atoms, name, values)


def _loader_for_trajectory(trajectory: _FakeTrajectory) -> TrajectoryLoader:
    """Build a loader instance backed by a fake trajectory."""

    loader = TrajectoryLoader.__new__(TrajectoryLoader)
    loader.load_universe = MagicMock(return_value=SimpleNamespace(trajectory=trajectory))
    return loader


class TestTrajectoryTimingMetadata:
    """Loader timing helpers preserve timestamps and reader position."""

    def test_get_first_frame_time_returns_ps_and_ns(self) -> None:
        """First-frame timestamps should be converted from MDAnalysis picoseconds."""

        trajectory = _FakeTrajectory([198_400.0, 198_800.0], current_frame=1)
        loader = _loader_for_trajectory(trajectory)

        assert loader.get_first_frame_time(1, unit="ps") == pytest.approx(198_400.0)
        assert loader.get_first_frame_time(1, unit="ns") == pytest.approx(198.4)

    @pytest.mark.parametrize(
        ("times_ps", "expose_time"),
        [([float("nan"), 400.0], True), ([0.0, 400.0], False)],
    )
    def test_get_first_frame_time_returns_none_for_unusable_time(
        self,
        times_ps: list[float],
        expose_time: bool,
    ) -> None:
        """Non-finite or missing trajectory times should return ``None``."""

        trajectory = _FakeTrajectory(times_ps, current_frame=1, expose_time=expose_time)
        loader = _loader_for_trajectory(trajectory)

        assert loader.get_first_frame_time(1) is None

    def test_get_first_frame_time_restores_current_frame(self) -> None:
        """Timestamp probing should restore the previous reader frame when feasible."""

        trajectory = _FakeTrajectory([198_400.0, 198_800.0, 199_200.0], current_frame=2)
        loader = _loader_for_trajectory(trajectory)

        assert loader.get_first_frame_time(1) == pytest.approx(198_400.0)
        assert trajectory.frame == 2

    def test_get_first_frame_time_prefers_raw_timestep_metadata(self) -> None:
        """Raw timestep times should win over normalized ChainReader times."""

        trajectory = _FakeTrajectory(
            [0.0, 400.0],
            raw_times_ps=[198_400.0, 198_800.0],
            current_frame=1,
        )
        loader = _loader_for_trajectory(trajectory)

        assert loader.get_first_frame_time(1) == pytest.approx(198_400.0)
        assert trajectory.frame == 1

    def test_get_frame_times_prefers_raw_timestep_metadata_and_restores_frame(self) -> None:
        """Frame-time extraction should use raw timestamps and restore reader position."""

        trajectory = _FakeTrajectory(
            [0.0, 400.0, 800.0],
            raw_times_ps=[198_400.0, 198_800.0, 199_200.0],
            current_frame=1,
        )
        loader = _loader_for_trajectory(trajectory)

        times = loader.get_frame_times(1, unit="ps")

        assert times.tolist() == pytest.approx([198_400.0, 198_800.0, 199_200.0])
        assert trajectory.frame == 1

    def test_get_timestep_restores_current_frame(self) -> None:
        """Timestep probing should restore the previous reader frame when feasible."""

        trajectory = _FakeTrajectory([198_400.0, 198_800.0, 199_200.0], current_frame=2)
        loader = _loader_for_trajectory(trajectory)

        assert loader.get_timestep(1, unit="ps") == pytest.approx(400.0)
        assert trajectory.frame == 2


class TestUniverseElementEnrichment:
    """Universe element enrichment handles GRO-like missing element metadata."""

    def test_existing_elements_are_preserved(self) -> None:
        """Existing element metadata should not be overwritten."""

        universe = _ElementFakeUniverse(
            _ElementFakeAtoms(names=["CA", "H1"], types=["C", "H"], elements=["C", "H"])
        )

        metadata = enrich_universe_elements(universe, topology_key="existing.gro")

        assert metadata["applied"] is False
        assert metadata["source"] == "existing"
        assert universe.added_topology_attrs == {}

    def test_safe_atom_types_add_elements(self) -> None:
        """Element-like GROMACS atom types should be copied as elements."""

        universe = _ElementFakeUniverse(
            _ElementFakeAtoms(names=["CA", "H1", "NA", "CL"], types=["C", "H", "NA", "CL"])
        )

        metadata = enrich_universe_elements(universe, topology_key="safe-types.gro")

        assert metadata == {
            "applied": True,
            "source": "types",
            "message": "elements inferred from atom types",
        }
        assert universe.atoms.elements == ["C", "H", "Na", "Cl"]

    def test_unsafe_atom_types_fall_back_to_names(self) -> None:
        """Force-field-like atom types should use conservative atom names."""

        universe = _ElementFakeUniverse(
            _ElementFakeAtoms(
                names=["CA", "1HA", "OW", "HW1"],
                types=["CT", "HC", "OW", "HW1"],
                resnames=["ALA", "ALA", "HOH", "HOH"],
            )
        )

        metadata = enrich_universe_elements(universe, topology_key="unsafe-types.gro")

        assert metadata["applied"] is True
        assert metadata["source"] == "names"
        assert universe.atoms.elements == ["C", "H", "O", "H"]

    @pytest.mark.parametrize(
        ("name", "resname", "expected"),
        [
            ("CA", "ALA", "C"),
            ("CA", "CAL", "Ca"),
            ("NA", "SOD", "Na"),
            ("CL", "CLA", "Cl"),
        ],
    )
    def test_ion_aliases_are_inferred_from_names(
        self,
        name: str,
        resname: str,
        expected: str,
    ) -> None:
        """Common ion aliases should distinguish ions from biomolecular names."""

        universe = _ElementFakeUniverse(
            _ElementFakeAtoms(names=[name], types=["XX"], resnames=[resname])
        )

        metadata = enrich_universe_elements(universe, topology_key=f"ion-{name}-{resname}.gro")

        assert metadata["applied"] is True
        assert metadata["source"] == "names"
        assert universe.atoms.elements == [expected]

    @pytest.mark.parametrize(
        ("name", "resname"),
        [("NA", "CLA"), ("CAX", "CAL")],
    )
    def test_ion_residue_conflicts_skip_enrichment(self, name: str, resname: str) -> None:
        """Recognized ion residues should reject conflicting atom-name aliases."""

        universe = _ElementFakeUniverse(
            _ElementFakeAtoms(names=[name], types=["XX"], resnames=[resname])
        )

        metadata = enrich_universe_elements(
            universe, topology_key=f"ion-conflict-{name}-{resname}.gro"
        )

        assert metadata["applied"] is False
        assert metadata["source"] is None
        assert not hasattr(universe.atoms, "elements")

    def test_unguessable_names_skip_enrichment(self) -> None:
        """No partial elements should be added when any atom is ambiguous."""

        universe = _ElementFakeUniverse(
            _ElementFakeAtoms(names=["CA", "BB"], types=["CT", "XX"], resnames=["ALA", "UNK"])
        )

        metadata = enrich_universe_elements(universe, topology_key="unguessable.gro")

        assert metadata["applied"] is False
        assert metadata["source"] is None
        assert not hasattr(universe.atoms, "elements")

    @pytest.mark.skipif(
        not Path("/home/joelaforet/Shirts-Lab-Linux/polyzymd_debugging/fn-d10.gro").exists(),
        reason="external GRO smoke-test topology is not available",
    )
    def test_external_gro_smoke_infers_elements(self) -> None:
        """External beta GRO topology should gain element metadata when present."""

        import MDAnalysis as mda

        gro_path = Path("/home/joelaforet/Shirts-Lab-Linux/polyzymd_debugging/fn-d10.gro")
        universe = mda.Universe(str(gro_path))

        metadata = enrich_universe_elements(universe, topology_key=gro_path)

        assert metadata["applied"] is True
        assert len(universe.atoms.elements) == len(universe.atoms)

    @pytest.mark.skipif(
        not Path(
            "/home/joelaforet/Shirts-Lab-Linux/polyzymd_debugging/fn-d10_OEG_n25.gro"
        ).exists(),
        reason="external polymer GRO smoke-test topology is not available",
    )
    def test_external_polymer_gro_enables_element_h_selection(self) -> None:
        """Element enrichment should preserve GRO hydrogen-selection UX for polymers."""

        import MDAnalysis as mda

        gro_path = Path("/home/joelaforet/Shirts-Lab-Linux/polyzymd_debugging/fn-d10_OEG_n25.gro")
        union_selection = "(protein) or (resid 95:219 and not resname HOH)"
        name_hydrogen_selection = f"({union_selection}) and (name H* or name [123]H*)"
        element_hydrogen_selection = f"({union_selection}) and element H"

        raw_universe = mda.Universe(str(gro_path))
        polymer_atoms = raw_universe.select_atoms("resid 95:219 and not resname HOH")
        name_hydrogens = raw_universe.select_atoms(name_hydrogen_selection)

        assert len(polymer_atoms) > 0
        assert len(name_hydrogens) == 4400
        if hasattr(raw_universe.atoms, "elements"):
            with pytest.raises(AttributeError):
                raw_universe.select_atoms(element_hydrogen_selection)
        else:
            with pytest.raises(AttributeError):
                raw_universe.select_atoms("element H")

        enriched_universe = mda.Universe(str(gro_path))
        metadata = enrich_universe_elements(enriched_universe, topology_key=gro_path)
        element_hydrogens = enriched_universe.select_atoms(element_hydrogen_selection)

        assert metadata["applied"] is True
        assert len(enriched_universe.atoms.elements) == len(enriched_universe.atoms)
        assert len(element_hydrogens) == len(name_hydrogens)


# ---------------------------------------------------------------------------
# Lazy engine creation
# ---------------------------------------------------------------------------


class TestLazyEngine:
    """Engine should not be created until the first resolution call."""

    def test_engine_not_created_at_init(self, tmp_path):
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)
        assert loader._engine is None

    def test_engine_created_on_first_find_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        loader.find_topology(run_dir)
        assert loader._engine is not None

    def test_engine_cached_across_calls(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        loader.find_topology(run_dir)
        first_engine = loader._engine
        loader.find_topology(run_dir)
        assert loader._engine is first_engine


# ---------------------------------------------------------------------------
# OpenMM layouts
# ---------------------------------------------------------------------------


class TestOpenMMDaisyChain:
    """Loader resolves OpenMM daisy-chain segmented layouts."""

    def test_finds_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=2)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file == run_dir / "solvated_system.pdb"

    def test_finds_all_segments_in_order(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=3)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.n_segments == 3
        for i, f in enumerate(info.trajectory_files):
            assert f.name == f"production_{i}_trajectory.dcd"

    def test_rejects_zero_byte_segments(self, tmp_path):
        """Loader should surface empty canonical OpenMM trajectory segments."""
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=3)
        empty_segment = run_dir / "production_3"
        empty_segment.mkdir()
        (empty_segment / "production_3_trajectory.dcd").write_bytes(b"")
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="Empty OpenMM trajectory segment"):
            loader.get_trajectory_info(replicate=1)

    def test_replicate_routing(self, tmp_path):
        for rep in (1, 3):
            run_dir = tmp_path / f"run_{rep}"
            _create_openmm_daisy_chain(run_dir, n_segments=1)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info1 = loader.get_trajectory_info(replicate=1)
        info3 = loader.get_trajectory_info(replicate=3)
        assert info1.working_directory != info3.working_directory
        assert info1.replicate == 1
        assert info3.replicate == 3


class TestOpenMMNoncanonical:
    """Loader resolves non-canonical single-directory OpenMM layout."""

    def test_finds_noncanonical_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_noncanonical(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file.name == "production_topology.pdb"

    def test_finds_noncanonical_trajectory(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_noncanonical(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert len(info.trajectory_files) == 1
        assert info.trajectory_files[0].name == "production_trajectory.dcd"


# ---------------------------------------------------------------------------
# GROMACS layouts
# ---------------------------------------------------------------------------


class TestGromacsFlat:
    """Loader resolves flat GROMACS production layout."""

    def test_finds_pdb_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file.suffix == ".pdb"
        assert info.topology_format == "pdb"

    def test_finds_prod_xtc(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert len(info.trajectory_files) == 1
        assert info.trajectory_files[0].name == "prod.xtc"
        assert info.trajectory_format == "xtc"

    def test_n_segments_is_one(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.n_segments == 1

    def test_records_gro_warning_in_metadata(self, tmp_path):
        """GRO topology warnings should be available to provenance wrappers."""
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir, use_gro=True)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)

        assert info.topology_format == "gro"
        assert any("GRO topology" in warning for warning in info.warnings)


# ---------------------------------------------------------------------------
# GRO topology warning
# ---------------------------------------------------------------------------


class TestGROWarning:
    """GRO topology triggers a chain-ID warning exactly once."""

    def test_warns_on_gro_topology(self, tmp_path, caplog):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir, use_gro=True)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.shared.loader"):
            loader.find_topology(run_dir)

        assert any(
            "GRO topology" in record.message and record.name == "polyzymd.analyses.shared.loader"
            for record in caplog.records
        )

    def test_warns_only_once_for_same_path(self, tmp_path, caplog):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir, use_gro=True)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.shared.loader"):
            loader.find_topology(run_dir)
            loader.find_topology(run_dir)

        gro_warnings = [
            record.message
            for record in caplog.records
            if "GRO topology" in record.message and record.name == "polyzymd.analyses.shared.loader"
        ]
        assert len(gro_warnings) == 1

    def test_no_warning_for_pdb_topology(self, tmp_path, caplog):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir, use_gro=False)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.shared.loader"):
            loader.find_topology(run_dir)

        gro_warnings = [m for m in caplog.messages if "GRO topology" in m]
        assert len(gro_warnings) == 0


# ---------------------------------------------------------------------------
# engine_override parameter
# ---------------------------------------------------------------------------


class TestEngineOverride:
    """engine_override forces a specific engine regardless of config.engine."""

    def test_override_gromacs_on_openmm_config(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        # Config says openmm, but override says gromacs
        config = _make_openmm_config(tmp_path, engine="openmm")
        loader = TrajectoryLoader(config, engine_override="gromacs")

        info = loader.get_trajectory_info(replicate=1)
        assert info.trajectory_files[0].name == "prod.xtc"

    def test_override_openmm_on_gromacs_config(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=2)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config, engine_override="openmm")

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file.name == "solvated_system.pdb"
        assert len(info.trajectory_files) == 2

    def test_override_openmm_succeeds_when_config_engine_missing(self, tmp_path):
        """OpenMM override should bypass a missing config.engine value."""
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)
        config = SimpleNamespace(
            get_working_directory=lambda rep: tmp_path / f"run_{rep}",
            output=SimpleNamespace(effective_scratch_directory=tmp_path),
        )
        loader = TrajectoryLoader(config, engine_override="openmm")

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file == run_dir / "solvated_system.pdb"


# ---------------------------------------------------------------------------
# find_topology backward compatibility
# ---------------------------------------------------------------------------


class TestFindTopologyCompat:
    """find_topology(working_dir) works for direct plugin callers."""

    def test_arbitrary_working_dir(self, tmp_path):
        """Plugin passes an explicit working_dir, not from config."""
        arbitrary_dir = tmp_path / "some_other_dir"
        arbitrary_dir.mkdir()
        (arbitrary_dir / "solvated_system.pdb").write_text("ATOM")

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        result = loader.find_topology(arbitrary_dir)
        assert result == arbitrary_dir / "solvated_system.pdb"

    def test_raises_when_no_topology(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="No topology file found"):
            loader.find_topology(empty_dir)


# ---------------------------------------------------------------------------
# _infer_replicate edge cases
# ---------------------------------------------------------------------------


class TestInferReplicate:
    """Best-effort replicate parsing from directory names."""

    def test_parses_run_1(self, tmp_path):
        result = TrajectoryLoader._infer_replicate(tmp_path / "run_1")
        assert result == 1

    def test_parses_run_42(self, tmp_path):
        result = TrajectoryLoader._infer_replicate(tmp_path / "run_42")
        assert result == 42

    def test_fallback_for_unknown_dir_name(self, tmp_path):
        result = TrajectoryLoader._infer_replicate(tmp_path / "my_simulation")
        assert result == 1

    def test_handles_non_path_gracefully(self):
        """MagicMock or other non-Path objects fall back to 1."""
        result = TrajectoryLoader._infer_replicate(MagicMock())
        assert result == 1


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------


class TestErrorPaths:
    """FileNotFoundError raised for missing files or directories."""

    def test_missing_working_directory(self, tmp_path):
        config = _make_openmm_config(tmp_path)
        (tmp_path / "run_1").mkdir()
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="Working directory not found") as exc_info:
            loader.get_trajectory_info(replicate=99)
        message = str(exc_info.value)
        assert "Replicate: 99" in message
        assert f"Expected working directory: {tmp_path / 'run_99'}" in message
        assert f"Scratch directory: {tmp_path}" in message
        assert "Available replicates: 1" in message
        assert "Action:" in message

    def test_no_topology_in_get_trajectory_info(self, tmp_path):
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        # Trajectory but no topology
        prod_dir = run_dir / "production_0"
        prod_dir.mkdir()
        (prod_dir / "production_0_trajectory.dcd").write_bytes(b"\x00")

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="No topology file found"):
            loader.get_trajectory_info(replicate=1)

    def test_no_trajectories_in_get_trajectory_info(self, tmp_path):
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "solvated_system.pdb").write_text("ATOM")
        # Topology but no trajectories

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(
            FileNotFoundError,
            match="No production trajectory files found",
        ) as exc_info:
            loader.get_trajectory_info(replicate=1)
        message = str(exc_info.value)
        assert "Replicate: 1" in message
        assert f"Expected working directory: {run_dir}" in message
        assert "Action:" in message

    def test_find_trajectories_raises_on_empty_dir(self, tmp_path):
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="No production trajectory files found"):
            loader._find_trajectories(run_dir)


# ---------------------------------------------------------------------------
# Engine resolution errors
# ---------------------------------------------------------------------------


class TestEngineResolutionErrors:
    """Invalid engine settings propagate instead of falling back."""

    def test_unknown_engine_error_propagates(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)

        config = MagicMock()
        config.engine = "nonexistent_engine"
        config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
        config.output.effective_scratch_directory = tmp_path

        loader = TrajectoryLoader(config)

        with pytest.raises(ValueError, match="Unknown engine"):
            loader.find_topology(run_dir)

    def test_missing_engine_error_propagates(self, tmp_path):
        """Configs without an engine value should be rejected by the loader."""
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)

        config = SimpleNamespace(
            get_working_directory=lambda rep: tmp_path / f"run_{rep}",
            output=SimpleNamespace(effective_scratch_directory=tmp_path),
        )

        loader = TrajectoryLoader(config)

        with pytest.raises(ValueError, match="non-empty string engine"):
            loader.find_topology(run_dir)

    @pytest.mark.parametrize("engine_value", [None, "", "   ", 123])
    def test_invalid_engine_values_propagate(self, tmp_path, engine_value: object):
        """Null, empty, and non-string config.engine values should be rejected."""
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)
        config = MagicMock()
        config.engine = engine_value
        config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
        config.output.effective_scratch_directory = tmp_path

        loader = TrajectoryLoader(config)

        with pytest.raises(ValueError, match="non-empty string engine"):
            loader.find_topology(run_dir)


# ---------------------------------------------------------------------------
# TrajectoryInfo validation
# ---------------------------------------------------------------------------


class TestTrajectoryInfoValidation:
    """TrajectoryInfo.validate() checks file existence."""

    def test_validate_passes_with_real_files(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        info.validate()  # Should not raise

    def test_validate_fails_missing_topology(self, tmp_path):
        from polyzymd.analyses.shared.loader import TrajectoryInfo

        info = TrajectoryInfo(
            topology_file=tmp_path / "nonexistent.pdb",
            trajectory_files=[],
        )
        with pytest.raises(FileNotFoundError, match="Topology not found"):
            info.validate()
