"""GROMACS smoke test for the hydrogen-bonds plugin."""

from __future__ import annotations

import sys
import types
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.hydrogen_bonds import HydrogenBondsAnalysis, HydrogenBondSettings
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class _Atom:
    def __init__(self, chain_id: str, resid: int, resname: str, resindex: int) -> None:
        self.segid = chain_id
        self.chainID = chain_id
        self.resid = resid
        self.resname = resname
        self.resindex = resindex


class _AtomGroup:
    def __init__(self, indices: list[int]) -> None:
        self.indices = np.asarray(indices, dtype=int)

    def __len__(self) -> int:
        return int(self.indices.size)

    def __or__(self, other: "_AtomGroup") -> "_AtomGroup":
        merged = sorted(set(self.indices.tolist()) | set(other.indices.tolist()))
        return _AtomGroup(merged)

    def __and__(self, other: "_AtomGroup") -> "_AtomGroup":
        overlap = sorted(set(self.indices.tolist()) & set(other.indices.tolist()))
        return _AtomGroup(overlap)


class TestHydrogenBondsGromacsSmoke:
    """Smoke test for hydrogen bonds with GROMACS trajectory layout."""

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run hydrogen-bonds compute on a real GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = HydrogenBondsAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = HydrogenBondSettings()

        universe = make_mock_universe(n_frames=5, n_atoms=2)
        universe.atoms = {
            0: _Atom(chain_id="A", resid=10, resname="SER", resindex=0),
            1: _Atom(chain_id="C", resid=100, resname="OEG", resindex=1),
        }
        universe.select_atoms = MagicMock(
            side_effect=lambda selection, updating=True: {
                "chainid A": _AtomGroup([0]),
                "chainid C": _AtomGroup([1]),
            }[selection]
        )

        hbonds_array = np.array([[0, 0, 10, 1, 2.8, 160.0]], dtype=float)
        hbond_cls = MagicMock()
        hbond_instance = MagicMock()
        hbond_instance.results = types.SimpleNamespace(hbonds=hbonds_array)
        hbond_instance.run.return_value = None
        hbond_cls.return_value = hbond_instance

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=universe)
        mock_modules = {
            "MDAnalysis": fake_mda,
            "MDAnalysis.analysis": types.ModuleType("MDAnalysis.analysis"),
            "MDAnalysis.analysis.hydrogenbonds": types.ModuleType(
                "MDAnalysis.analysis.hydrogenbonds"
            ),
            "MDAnalysis.analysis.hydrogenbonds.hbond_analysis": types.ModuleType(
                "MDAnalysis.analysis.hydrogenbonds.hbond_analysis"
            ),
        }
        mock_modules["MDAnalysis.analysis.hydrogenbonds.hbond_analysis"].HydrogenBondAnalysis = (
            hbond_cls
        )

        original_resolve = GromacsEngine.resolve_trajectory_layout
        with (
            patch.dict(sys.modules, mock_modules),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
            patch("polyzymd.analyses.hydrogen_bonds.compute_config_hash", return_value="smoke123"),
        ):
            ctx = ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=config,
                output_dir=tmp_path / "analysis" / "run_1",
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            result = analysis.run_replicate(ctx, replicate=1)

        assert resolve_spy.call_count >= 1
        assert result.replicate == 1
        assert hbond_cls.call_args.kwargs["donors_sel"] == "(chainid A) or (chainid C)"
        assert str(tmp_path / "run_1" / "gromacs" / "prod.xtc") in str(fake_mda.Universe.call_args)
