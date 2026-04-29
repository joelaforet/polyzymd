"""GROMACS smoke test for the secondary-structure plugin."""

from __future__ import annotations

import sys
import types
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.secondary_structure import (
    SecondaryStructureAnalysis,
    SecondaryStructureSettings,
)
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestSecondaryStructureGromacsSmoke:
    """Smoke test for secondary structure with GROMACS layout resolution."""

    def test_compute_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run secondary-structure compute on a GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = SecondaryStructureAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = SecondaryStructureSettings(chain_id="A")

        mock_mdtraj_traj = MagicMock()
        mock_mdtraj_traj.n_frames = 20000
        mock_mdtraj_traj.n_atoms = 100
        mock_mdtraj_traj.n_residues = 5
        mock_mdtraj_traj.topology = MagicMock()
        mock_mdtraj_traj.topology.select.return_value = np.arange(20)
        mock_mdtraj_traj.topology.chains = [MagicMock(index=0)]
        mock_mdtraj_traj.topology.residues = [
            types.SimpleNamespace(resSeq=1, name="ALA"),
            types.SimpleNamespace(resSeq=2, name="GLY"),
            types.SimpleNamespace(resSeq=3, name="SER"),
            types.SimpleNamespace(resSeq=4, name="LEU"),
            types.SimpleNamespace(resSeq=5, name="VAL"),
        ]

        protein_traj = MagicMock()
        protein_traj.n_frames = 20000
        protein_traj.n_atoms = 20
        protein_traj.n_residues = 5
        protein_traj.topology = mock_mdtraj_traj.topology
        protein_traj.__getitem__.return_value = protein_traj
        mock_mdtraj_traj.atom_slice.return_value = protein_traj

        dssp_raw = np.array([list("HECHC")] * protein_traj.n_frames)
        ss_matrix = np.array([[1, 2, 0, 1, 0]] * protein_traj.n_frames, dtype=np.int8)

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=20000, n_atoms=50))

        original_resolve = GromacsEngine.resolve_trajectory_layout
        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
            patch(
                "polyzymd.analyses.secondary_structure.compute_config_hash", return_value="smoke123"
            ),
            patch("mdtraj.load", return_value=mock_mdtraj_traj) as mdtraj_load,
            patch("mdtraj.compute_dssp", return_value=dssp_raw),
            patch(
                "polyzymd.analyses.secondary_structure._encode_dssp_matrix", return_value=ss_matrix
            ),
            patch(
                "polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.3.0-test"
            ),
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

        assert resolve_spy.call_count >= 2
        assert result.replicate == 1
        assert result.selection_string == "chainid 0 (chain A)"
        assert result.overall_helix_fraction == 0.4
        assert result.overall_strand_fraction == 0.2
        assert str(tmp_path / "run_1" / "gromacs" / "prod.xtc") in str(mdtraj_load.call_args)
