"""Compare PolyzyMD's per-residue RMSF and RMS deviation with GROMACS ``gmx rmsf``.

For each condition and reference mode, PolyzyMD builds the reference with
``pz.reference`` and computes ``rmsf``, ``residue_rmsd`` and ``offset`` with
``polyzymd.analyses.functions.rms_decomposition``. The same production frames
of the same atoms, not superposed, are written to a TRR file, and the
reference positions to a GROMOS96 file, which stores 1e-4 nm. ``gmx rmsf -res``
then fits every frame to that reference and writes each residue's fluctuation
about its mean position (``-o``) and its deviation from the reference
(``-od``). The GROMACS offset is sqrt(deviation^2 - RMSF^2). ``gmx rmsf``
writes its .xvg values to 4 decimals in nm, so differences up to 5e-4 Å are
output rounding.

Use a selection with one atom per residue, such as Cα atoms: GROMACS then
fits with equal masses, as PolyzyMD does.

Run in the analysis pixi environment with ``gmx`` on the PATH::

    python scripts/benchmarks/rmsf_gromacs_parity.py \
        --config noPoly=path/to/noPoly/config.yaml --eq noPoly=200ns \
        --config sbma=path/to/sbma/config.yaml --eq sbma=10ns \
        --selection "protein and name CA and resid 4:174" \
        --reference-file structures/1ISP.pdb --frame 50 --output-dir gmx_parity

It prints one line per condition and mode and writes ``parity.json``.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import MDAnalysis as mda
import numpy as np

import polyzymd as pz
from polyzymd.analyses import functions
from polyzymd.analyses.reference import build_reference


def write_g96(path: Path, atoms, positions_angstrom: np.ndarray) -> None:
    """Write ``atoms`` at ``positions_angstrom`` as a GROMOS96 file, in nm."""
    lines = ["TITLE", "reference", "END", "POSITION"]
    for index, (atom, xyz) in enumerate(zip(atoms, positions_angstrom / 10.0), start=1):
        lines.append(
            f"{atom.resid:5d} {atom.resname:<5s} {atom.name:<5s}{index:7d}"
            f"{xyz[0]:15.9f}{xyz[1]:15.9f}{xyz[2]:15.9f}"
        )
    lines += ["END", "BOX", f"{20.0:15.9f}{20.0:15.9f}{20.0:15.9f}", "END"]
    path.write_text("\n".join(lines) + "\n")


def read_xvg(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return the residue IDs and values of a ``gmx rmsf`` .xvg file, values in Å."""
    rows = [line.split() for line in path.read_text().splitlines() if line and line[0] not in "#@"]
    data = np.array(rows, float)
    return data[:, 0].astype(int), data[:, 1] * 10.0


def compare(ours: np.ndarray, theirs: np.ndarray) -> dict[str, float]:
    """Pearson r, mean and maximum absolute difference of two per-residue profiles."""
    difference = np.abs(ours - theirs)
    spread = np.std(ours) > 0 and np.std(theirs) > 0
    return {
        "pearson_r": float(np.corrcoef(ours, theirs)[0, 1]) if spread else float("nan"),
        "mean_abs": float(difference.mean()),
        "max_abs": float(difference.max()),
        "mean_value": float(ours.mean()),
    }


def pairs(values: list[str], flag: str) -> dict[str, str]:
    """Parse ``name=value`` arguments."""
    parsed = dict(value.split("=", 1) for value in values)
    if len(parsed) != len(values):
        raise SystemExit(f"{flag}: repeated names in {values}")
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", action="append", required=True, help="name=config.yaml")
    parser.add_argument("--eq", action="append", required=True, help="name=equilibration")
    parser.add_argument("--selection", default="protein and name CA")
    parser.add_argument("--reference-file", required=True, help="Structure for external mode")
    parser.add_argument("--frame", type=int, default=1, help="Production frame for frame mode")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=Path("gmx_parity"))
    parser.add_argument("--gmx", default="gmx")
    args = parser.parse_args()

    configs, windows = pairs(args.config, "--config"), pairs(args.eq, "--eq")
    modes = {
        "external": {"file": args.reference_file},
        "average": {},
        "centroid": {},
        "frame": {"frame": args.frame},
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, config in configs.items():
        study = pz.Study.from_configs(
            {name: config}, equilibration=windows[name], replicates=[args.replicate]
        )
        replicate = study[name].replicates[0]
        universe, frames = replicate.universe(), replicate.frames
        atoms = universe.select_atoms(args.selection)
        for mode, extra in modes.items():
            reference = pz.reference(mode, args.selection, alignment=args.selection, **extra)
            ref_atoms, _ = build_reference(reference, universe, frames)
            deviation, fluctuation, offset = functions.rms_decomposition(
                atoms, atoms, ref_atoms, frames
            )[:3]

            work = args.output_dir / f"{name}_{mode}"
            work.mkdir(exist_ok=True)
            write_g96(work / "reference.g96", atoms, ref_atoms.positions.astype(float))
            with mda.Writer(str(work / "frames.trr"), atoms.n_atoms) as writer:
                for _ in universe.trajectory[frames]:
                    writer.write(atoms)
            command = [
                args.gmx,
                "rmsf",
                "-s",
                "reference.g96",
                "-f",
                "frames.trr",
                "-o",
                "rmsf.xvg",
                "-od",
                "deviation.xvg",
                "-res",
                "-quiet",
            ]
            run = subprocess.run(command, cwd=work, input="0\n", capture_output=True, text=True)
            if run.returncode:
                raise SystemExit(f"gmx rmsf failed for {name} {mode}:\n{run.stderr[-2000:]}")
            residues, gmx_fluctuation = read_xvg(work / "rmsf.xvg")
            _, gmx_deviation = read_xvg(work / "deviation.xvg")
            if list(residues) != list(atoms.residues.resids):
                raise SystemExit(f"{name} {mode}: GROMACS residues differ from the selection")
            gmx_offset = np.sqrt(np.clip(gmx_deviation**2 - gmx_fluctuation**2, 0.0, None))

            row = {"n_frames": int(len(frames)), "n_residues": int(len(residues))}
            row["rmsf"] = compare(fluctuation, gmx_fluctuation)
            row["residue_rmsd"] = compare(deviation, gmx_deviation)
            row["offset"] = compare(offset, gmx_offset)
            results[f"{name} {mode}"] = row
            print(
                f"{name} {mode}: frames {row['n_frames']}, residues {row['n_residues']} | "
                + " | ".join(
                    f"{part} r {row[part]['pearson_r']:.8f} mean {row[part]['mean_abs']:.1e} "
                    f"max {row[part]['max_abs']:.1e} A"
                    for part in ("rmsf", "residue_rmsd", "offset")
                ),
                flush=True,
            )
    version = subprocess.run([args.gmx, "--version"], capture_output=True, text=True).stdout
    gmx_version = next(
        (
            line.split(":", 1)[1].strip()
            for line in version.splitlines()
            if "GROMACS version" in line
        ),
        "?",
    )
    results["versions"] = {
        "gromacs": gmx_version,
        "MDAnalysis": mda.__version__,
        "polyzymd": pz.__version__,
    }
    (args.output_dir / "parity.json").write_text(json.dumps(results, indent=1))
    print(f"GROMACS {gmx_version}, MDAnalysis {mda.__version__}, PolyzyMD {pz.__version__}")


if __name__ == "__main__":
    main()
