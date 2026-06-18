#!/usr/bin/env python
"""Subprocess script for polymer generation using polyzymd build env.

This script runs under the polyzymd build environment (Python 3.11 + mbuild)
to generate 3D polymer structures with correct monomer placement.

It fixes the polymerist block identifier mapping bug by passing an explicit
``sequence_map`` to ``build_linear_polymer()``, bypassing the default
alphabetical-sort-based mapping that disagrees with reverse-sorted
MonomerGroup keys.

Usage (from conjugation env, via subprocess)::

    /path/to/polyzymd/build/python generate_polymer.py \\
        --cache-dir .polymer_cache_conjugation \\
        --monogrp-json .polymer_cache_conjugation/NHS-SBMA_monomer_group.json \\
        --sequence BAB \\
        --monomer-names '{"A": "NHS", "B": "SBMA"}' \\
        --residue-names '{"NHS": "NH", "SBMA": "SB"}'

Outputs:
    - <cache_dir>/<prefix>_seq=<seq>_<n>-mer.pdb
    - <cache_dir>/<prefix>_seq=<seq>_<n>-mer.sdf  (uncharged)

Prints the PDB path to stdout on success.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("generate_polymer")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate polymer 3D structure")
    parser.add_argument("--cache-dir", required=True, help="Cache directory for output")
    parser.add_argument("--monogrp-json", required=True, help="Path to MonomerGroup JSON")
    parser.add_argument("--sequence", required=True, help="Full sequence string e.g. BAB")
    parser.add_argument(
        "--monomer-names",
        required=True,
        help='JSON dict mapping sequence labels to monomer names, e.g. \'{"A": "NHS", "B": "SBMA"}\'',
    )
    parser.add_argument(
        "--residue-names",
        required=False,
        default="{}",
        help='JSON dict mapping monomer names to 2-char residue prefixes, e.g. \'{"NHS": "NH", "SBMA": "SB"}\'',
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=10,
        help="Max retries for ring-piercing avoidance",
    )
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    monomer_names: dict[str, str] = json.loads(args.monomer_names)
    residue_names: dict[str, str] = json.loads(args.residue_names)
    sequence: str = args.sequence

    # Load MonomerGroup from cached JSON
    from polymerist.polymers.monomers.repr import MonomerGroup

    monogrp = MonomerGroup.from_file(args.monogrp_json)
    logger.info("Loaded MonomerGroup with keys: %s", list(monogrp.monomers.keys()))

    # Parse sequence: head + middle + tail
    if len(sequence) < 2:
        raise ValueError(f"Sequence must be at least 2 characters, got '{sequence}'")

    head_label = sequence[0]
    tail_label = sequence[-1]
    middle_sequence = sequence[1:-1] if len(sequence) > 2 else ""

    head_monomer = monomer_names[head_label]
    tail_monomer = monomer_names[tail_label]

    # Helper to find the 1-site fragment name for a monomer
    def get_terminal_fragment_name(monomer_name: str) -> str:
        for key in monogrp.monomers:
            if key.startswith(monomer_name) and "_1-site" in key:
                return key
        raise ValueError(f"No 1-site fragment found for monomer '{monomer_name}'")

    # Helper to find the 2-site fragment name for a monomer
    def get_middle_fragment_name(monomer_name: str) -> str:
        for key in monogrp.monomers:
            if key.startswith(monomer_name) and "_2-site" in key:
                return key
        raise ValueError(f"No 2-site fragment found for monomer '{monomer_name}'")

    # Set terminal orientations
    monogrp_local = MonomerGroup(monomers=monogrp.monomers)
    monogrp_local.term_orient = {
        "head": get_terminal_fragment_name(head_monomer),
        "tail": get_terminal_fragment_name(tail_monomer),
    }
    logger.info(
        "Terminal orientations: head=%s, tail=%s",
        monogrp_local.term_orient["head"],
        monogrp_local.term_orient["tail"],
    )

    # Build explicit sequence_map: map each unique middle character to the
    # correct 2-site fragment name based on monomer_names dict
    if middle_sequence:
        unique_middle_chars = []
        for ch in middle_sequence:
            if ch not in unique_middle_chars:
                unique_middle_chars.append(ch)

        sequence_map = {}
        for ch in unique_middle_chars:
            if ch not in monomer_names:
                raise ValueError(
                    f"Middle sequence character '{ch}' not found in monomer_names {monomer_names}"
                )
            fragment_name = get_middle_fragment_name(monomer_names[ch])
            sequence_map[ch] = fragment_name

        logger.info("Explicit sequence_map: %s", sequence_map)
    else:
        sequence_map = None
        logger.info("No middle sequence; only terminals will be built")

    # Import build_linear_polymer and utilities (all from polymerist)
    from polymerist.polymers.building import mbmol_to_rdmol
    from polymerist.polymers.building.linear import build_linear_polymer
    from polymerist.rdutils.rdcoords.piercing import summarize_ring_piercing

    # Build with retries
    chain = None
    for attempt in range(args.max_retries):
        logger.info("Building polymer attempt %d/%d", attempt + 1, args.max_retries)
        try:
            chain = build_linear_polymer(
                monomers=monogrp_local,
                n_monomers=len(sequence),
                sequence=middle_sequence if middle_sequence else "A",
                sequence_map=sequence_map,
                energy_minimize=True,
                allow_partial_sequences=True,
            )
        except Exception as exc:
            logger.warning("Build failed on attempt %d: %s", attempt + 1, exc)
            continue

        # Check ring-piercing
        poly_mol = mbmol_to_rdmol(chain)
        piercing = summarize_ring_piercing(poly_mol)
        if not piercing:
            logger.info("Polymer built successfully on attempt %d", attempt + 1)
            break
        else:
            logger.warning("Ring-piercing on attempt %d: %s", attempt + 1, piercing)
            chain = None
    else:
        if chain is None:
            print("FAILED: Could not build polymer after max retries", file=sys.stderr)
            sys.exit(1)

    # Build resname_map for PDB output
    resname_map: dict[str, str] = {}
    for label, monomer_name in monomer_names.items():
        if monomer_name in residue_names:
            base = residue_names[monomer_name]
        else:
            base = monomer_name[:2].upper()
        for frag_name in monogrp.monomers:
            if frag_name.startswith(monomer_name):
                if "_1-site" in frag_name:
                    resname_map[frag_name] = f"{base}1"
                elif "_2-site" in frag_name:
                    resname_map[frag_name] = f"{base}2"
    logger.info("Residue name map: %s", resname_map)

    # Build filename
    unique_labels = []
    for ch in sequence:
        if ch not in unique_labels:
            unique_labels.append(ch)
    monomers_used = [monomer_names[ch] for ch in unique_labels]
    prefix = "-".join(monomers_used)
    filename = f"{prefix}_seq={sequence}_{len(sequence)}-mer"

    pdb_path = cache_dir / f"{filename}.pdb"
    sdf_path = cache_dir / f"{filename}.sdf"

    # Save PDB using polymerist's utility
    from polymerist.polymers.building import mbmol_to_openmm_pdb

    mbmol_to_openmm_pdb(pdb_path, chain, resname_map=resname_map)
    logger.info("Saved PDB: %s", pdb_path)

    # Save uncharged SDF via OpenFF topology
    try:
        from openff.toolkit import Topology as OFFTopology

        from polyzymd.utils.topology import topology_to_sdf

        off_top = OFFTopology.from_pdb(str(pdb_path), _custom_substructures=monogrp.monomers)
        # Partition into residues
        from polymerist.mdtools.openfftools.partition import partition

        was_partitioned = partition(off_top)
        if not was_partitioned:
            logger.warning("Partitioning failed, saving SDF without residue info")

        # Truncate extended residue names
        for mol in off_top.molecules:
            for atom in mol.atoms:
                if "residue_name" in atom.metadata:
                    atom.metadata["extended_name"] = atom.metadata["residue_name"]
                    atom.metadata["residue_name"] = atom.metadata["residue_name"][:3]

        topology_to_sdf(sdf_path, off_top)
        logger.info("Saved uncharged SDF: %s", sdf_path)
    except Exception as exc:
        logger.warning("SDF generation failed (non-fatal): %s", exc)

    # Print the PDB path to stdout for the caller
    print(str(pdb_path))


if __name__ == "__main__":
    main()
