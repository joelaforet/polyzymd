"""Check PolyzyMD's SASA against an independent Shrake-Rupley calculation, frame by frame.

For each of the first production frames of one replicate, the script computes
the SASA of ``--target`` in ``--context`` three ways:

- ``polyzymd.analyses.functions.sasa``, which gives MDTraj one frame per call;
- an independent NumPy/SciPy Shrake-Rupley calculation with the same atomic
  radii (MDTraj's element table plus the probe) and the same golden-section
  spiral of sphere points as MDTraj;
- ``mdtraj.shrake_rupley`` given the same frame ``--copies`` times in one call,
  which shows whether a frame's area depends on the frames computed before it
  in the same call.

Run in the analysis pixi environment::

    python scripts/benchmarks/sasa_mdtraj_check.py --config path/to/config.yaml \
        --eq 10ns --target protein --context "protein or resname SBM EGM" --frames 3

It prints one line per frame and writes ``sasa_check.json`` to ``--output-dir``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
from mdtraj.geometry.sasa import _ATOMIC_RADII
from scipy.spatial import cKDTree

import polyzymd as pz
from polyzymd.analyses import functions
from polyzymd.analyses.functions import _sasa_topology


def spiral(n_points: int) -> np.ndarray:
    """Return MDTraj's golden-section spiral of ``n_points`` unit vectors."""
    increment, offset = np.pi * (3.0 - np.sqrt(5.0)), 2.0 / n_points
    k = np.arange(n_points)
    y = k * offset - 1.0 + offset / 2.0
    r = np.sqrt(1.0 - y * y)
    phi = k * increment
    return np.stack([np.cos(phi) * r, y, np.sin(phi) * r], axis=1)


def shrake_rupley(xyz_nm: np.ndarray, radii_nm: np.ndarray, n_points: int) -> np.ndarray:
    """Return each atom's SASA in nm², counting the sphere points inside no other sphere."""
    points, tree, largest = spiral(n_points), cKDTree(xyz_nm), radii_nm.max()
    area = np.zeros(len(xyz_nm))
    for i, (centre, radius) in enumerate(zip(xyz_nm, radii_nm)):
        neighbours = np.array(
            [j for j in tree.query_ball_point(centre, radius + largest) if j != i], dtype=int
        )
        surface = centre + radius * points
        if neighbours.size:
            squared = ((surface[:, None, :] - xyz_nm[neighbours][None]) ** 2).sum(-1)
            free = ~(squared < radii_nm[neighbours][None] ** 2).any(axis=1)
        else:
            free = np.ones(n_points, bool)
        area[i] = 4.0 * np.pi * radius * radius * free.mean()
    return area


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    parser.add_argument("--eq", required=True, help="Equilibration window, for example 10ns")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--target", default="protein")
    parser.add_argument("--context", default=None, help="Defaults to the target")
    parser.add_argument("--frames", type=int, default=3, help="Production frames to check")
    parser.add_argument("--copies", type=int, default=3)
    parser.add_argument("--probe-radius-nm", type=float, default=functions.SASA_PROBE_RADIUS_NM)
    parser.add_argument("--n-sphere-points", type=int, default=functions.SASA_SPHERE_POINTS)
    parser.add_argument("--output-dir", type=Path, default=Path("sasa_check"))
    args = parser.parse_args()

    study = pz.Study.from_configs(
        {"c": args.config}, equilibration=args.eq, replicates=[args.replicate]
    )
    replicate = study["c"].replicates[0]
    universe = replicate.universe()
    target = universe.select_atoms(args.target)
    context = universe.select_atoms(args.context or args.target)
    rows = []
    for frame in replicate.frames[: args.frames]:
        universe.trajectory[frame]
        topology, where = _sasa_topology(target, context)
        xyz = (context.positions / 10.0).astype(np.float32)
        radii = (
            np.array([_ATOMIC_RADII[a.element.symbol] for a in topology.atoms])
            + args.probe_radius_nm
        )
        independent = (
            shrake_rupley(xyz.astype(np.float64), radii, args.n_sphere_points)[where].sum() * 100
        )
        polyzymd = functions.sasa(target, context, args.probe_radius_nm, args.n_sphere_points)
        repeated = (
            md.shrake_rupley(
                md.Trajectory(xyz=np.repeat(xyz[None], args.copies, axis=0), topology=topology),
                mode="atom",
                probe_radius=args.probe_radius_nm,
                n_sphere_points=args.n_sphere_points,
            )[:, where].sum(axis=1)
            * 100
        )
        rows.append(
            {
                "frame": int(frame),
                "independent": float(independent),
                "polyzymd": polyzymd,
                "mdtraj_repeated": [float(v) for v in repeated],
            }
        )
        print(
            f"frame {frame}: independent {independent:.3f}, polyzymd {polyzymd:.3f} A^2 "
            f"({100 * (polyzymd - independent) / independent:+.4f}%); MDTraj with {args.copies} "
            f"copies in one call {', '.join(f'{v:.3f}' for v in repeated)}",
            flush=True,
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    versions = {"mdtraj": md.__version__, "polyzymd": pz.__version__}
    payload = {"target": args.target, "context": args.context, "rows": rows, "versions": versions}
    (args.output_dir / "sasa_check.json").write_text(json.dumps(payload, indent=1))
    print(
        f"MDTraj {md.__version__}, PolyzyMD {pz.__version__}, {len(target)} target and {len(context)} context atoms"
    )


if __name__ == "__main__":
    main()
