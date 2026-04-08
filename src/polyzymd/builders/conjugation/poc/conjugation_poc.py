import marimo

__generated_with = "0.22.5"
app = marimo.App()


@app.cell
def _():
    """Cell 1: Configuration, imports, and shared utilities."""
    from dataclasses import dataclass, field
    from pathlib import Path
    import json
    import logging
    import subprocess

    import numpy as np
    from scipy.spatial import KDTree
    from scipy.spatial.transform import Rotation
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from openff.toolkit import ForceField, Molecule, Topology
    from openff.units import Quantity, unit

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("conjugation_poc")

    # ── Paths ──────────────────────────────────────────────────────────────
    PROTEIN_PDB = (
        "/home/joelaforet/Shirts-Lab-Linux/openff-pablo/SBMA_EGMA_50_50_water_363K/"
        "structures/NH3_terminal_His_proton_updated.pdb"
    )
    ASSEMBLED_PDB = Path(
        "/home/joelaforet/Shirts-Lab-Linux/2026-virtual-workshops/conjugate_assembled.pdb"
    )
    MINIMIZED_PDB = Path(
        "/home/joelaforet/Shirts-Lab-Linux/2026-virtual-workshops/conjugate_minimized.pdb"
    )
    POLYZYMD_PYTHON = (
        "/home/joelaforet/Shirts-Lab-Linux/polyzymd/.pixi/envs/build/bin/python"
    )
    GENERATE_SCRIPT = (
        "/home/joelaforet/Shirts-Lab-Linux/2026-virtual-workshops/generate_polymer.py"
    )
    POLYMER_CACHE = Path(
        "/home/joelaforet/Shirts-Lab-Linux/2026-virtual-workshops/.polymer_cache_conjugation"
    )
    MONOMER_GROUP_JSON = POLYMER_CACHE / "NHS-SBMA_monomer_group.json"

    # ── Conjugation parameters ─────────────────────────────────────────────
    ELIGIBLE_LYS_RESIDS = [23, 35, 44]
    N_CONJUGATED = 2
    POLYMER_SEQUENCE = "BAB"  # A=NHS, B=SBMA; first residue contains NHS
    MONOMER_NAMES = {"A": "NHS", "B": "SBMA"}
    RESIDUE_NAMES = {"NHS": "NH", "SBMA": "SB"}

    # ── Placement parameters ───────────────────────────────────────────────
    N_CONFORMERS = 50
    N_CONE_SAMPLES = 50
    N_SPIN_SAMPLES = 50
    COLLISION_DISTANCE = 1.5       # angstrom
    BOND_ANGLE_TOLERANCE = 25.0    # degrees
    TARGET_BOND_LENGTH = 1.33      # angstrom (amide C-N)
    TARGET_BOND_ANGLE = 120.0      # degrees (at NZ)
    RNG_SEED = 42

    # ── Force field / minimization ─────────────────────────────────────────
    PROTEIN_FF = "ff14sb_off_impropers_0.0.4.offxml"
    SMALL_MOL_FF = "openff-2.0.0.offxml"
    MINIMIZATION_MAX_ITER = 500
    LINKAGE_NEIGHBORHOOD_BONDS = 3

    # ── Data classes ───────────────────────────────────────────────────────
    @dataclass(frozen=True)
    class ReactiveSite:
        """Reactive lysine site on the protein.

        Attributes
        ----------
        resid : int
            PDB residue sequence number.
        chain_id : str
            PDB chain identifier.
        nz_index : int
            Atom index of NZ in the protein OpenFF Molecule.
        nz_position : np.ndarray
            Cartesian position of NZ (angstrom).
        ce_position : np.ndarray
            Cartesian position of CE (angstrom).
        hz_indices : tuple[int, ...]
            Atom indices of NZ-bound hydrogens in protein OpenFF Molecule.
        outward_normal : np.ndarray
            Unit vector from protein centroid toward NZ.
        """
        resid: int
        chain_id: str
        nz_index: int
        nz_position: np.ndarray
        ce_position: np.ndarray
        hz_indices: tuple[int, ...]
        outward_normal: np.ndarray

    @dataclass(frozen=True)
    class PolymerPrototype:
        """Multi-residue polymer oligomer loaded from polyzymd.

        Attributes
        ----------
        off_molecule : Molecule
            OpenFF Molecule with full connectivity and conformer.
        rdmol : Chem.Mol
            RDKit molecule (all-single-bond from PDB load).
        reactive_carbon_idx : int
            Index of NHS ester carbonyl carbon in the RDKit/OFF molecule.
        leaving_group_indices : tuple[int, ...]
            Indices of NHS leaving group atoms (ester O + succinimide ring).
        retained_atom_indices : tuple[int, ...]
            Indices of atoms kept after leaving group removal.
        product_charges : np.ndarray
            NAGL partial charges for the retained atoms (product state),
            indexed by retained_atom_indices order.
        """
        off_molecule: Molecule
        rdmol: object  # Chem.Mol
        reactive_carbon_idx: int
        leaving_group_indices: tuple[int, ...]
        retained_atom_indices: tuple[int, ...]
        product_charges: np.ndarray

    @dataclass(frozen=True)
    class PlacementResult:
        """Record of a successful geometric placement.

        Attributes
        ----------
        site : ReactiveSite
            Protein reactive site used.
        polymer : PolymerPrototype
            Polymer placed at this site.
        placed_coords : np.ndarray
            Full precursor coordinates after placement (angstrom).
        min_protein_distance : float
            Best minimum distance to occupied atoms (angstrom).
        """
        site: ReactiveSite
        polymer: PolymerPrototype
        placed_coords: np.ndarray
        min_protein_distance: float

    # ── PDB writing utilities ──────────────────────────────────────────────
    def write_pdb_line(
        serial: int,
        atom_name: str,
        residue_name: str,
        chain_id: str,
        residue_number: int,
        x: float,
        y: float,
        z: float,
        element: str,
        record: str = "HETATM",
    ) -> str:
        if len(atom_name) < 4:
            padded_name = f" {atom_name:<3}"
        else:
            padded_name = f"{atom_name:<4}"
        return (
            f"{record:<6}{serial:>5d} {padded_name} {residue_name:>3} {chain_id:1}"
            f"{residue_number:>4d}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.00:>6.2f}{0.00:>6.2f}          {element:>2}\n"
        )

    def write_conect_line(serial: int, bonded_serials: list[int]) -> str:
        lines: list[str] = []
        for start in range(0, len(bonded_serials), 4):
            chunk = bonded_serials[start : start + 4]
            fields = "".join(f"{entry:>5d}" for entry in chunk)
            lines.append(f"CONECT{serial:>5d}{fields}\n")
        return "".join(lines)

    return (
        ASSEMBLED_PDB,
        AllChem,
        BOND_ANGLE_TOLERANCE,
        COLLISION_DISTANCE,
        Chem,
        ELIGIBLE_LYS_RESIDS,
        ForceField,
        GENERATE_SCRIPT,
        KDTree,
        LINKAGE_NEIGHBORHOOD_BONDS,
        MINIMIZATION_MAX_ITER,
        MINIMIZED_PDB,
        MONOMER_GROUP_JSON,
        MONOMER_NAMES,
        Molecule,
        N_CONFORMERS,
        N_CONJUGATED,
        N_CONE_SAMPLES,
        N_SPIN_SAMPLES,
        POLYMER_CACHE,
        POLYMER_SEQUENCE,
        POLYZYMD_PYTHON,
        PROTEIN_FF,
        PROTEIN_PDB,
        Path,
        PlacementResult,
        PolymerPrototype,
        Quantity,
        RESIDUE_NAMES,
        RNG_SEED,
        ReactiveSite,
        Rotation,
        SMALL_MOL_FF,
        TARGET_BOND_ANGLE,
        TARGET_BOND_LENGTH,
        Topology,
        json,
        logger,
        np,
        subprocess,
        unit,
        write_conect_line,
        write_pdb_line,
    )


@app.cell
def _(
    ELIGIBLE_LYS_RESIDS,
    ForceField,
    PROTEIN_FF,
    PROTEIN_PDB,
    Quantity,
    ReactiveSite,
    Topology,
    logger,
    np,
    unit,
):
    """Cell 2: Load protein PDB, identify reactive lysine sites, extract ff14SB charges."""

    # Load protein as OpenFF Topology → single Molecule.
    # Topology.from_pdb uses the CCD residue library internally, same as pablo,
    # but we don't need pablo's crosslink/linking_bond machinery here.
    protein_top = Topology.from_pdb(PROTEIN_PDB)
    protein_mol = list(protein_top.molecules)[0]
    protein_coords = protein_mol.conformers[0].m_as("angstrom")
    centroid_vec = protein_coords.mean(axis=0)

    logger.info(
        "Loaded protein: %d atoms, %d bonds",
        protein_mol.n_atoms,
        protein_mol.n_bonds,
    )

    # ── Identify reactive lysine sites ─────────────────────────────────────
    reactive_sites: dict[tuple[str, int], ReactiveSite] = {}
    for _resid in ELIGIBLE_LYS_RESIDS:
        nz_atom = None
        for _atom in protein_mol.atoms:
            if (
                _atom.metadata.get("residue_number") == str(_resid)
                and _atom.metadata.get("residue_name") == "LYS"
                and _atom.name == "NZ"
            ):
                nz_atom = _atom
                break
        if nz_atom is None:
            raise ValueError(f"Eligible lysine {_resid} lacks NZ")

        ce_pos = None
        hz_idx_list: list[int] = []
        for _bond in nz_atom.bonds:
            other_idx = (
                _bond.atom1_index
                if _bond.atom2_index == nz_atom.molecule_atom_index
                else _bond.atom2_index
            )
            neighbor = protein_mol.atom(other_idx)
            if neighbor.name == "CE":
                ce_pos = protein_coords[other_idx]
            if neighbor.atomic_number == 1:
                hz_idx_list.append(other_idx)
        if ce_pos is None:
            raise ValueError(f"Lysine {_resid} lacks CE")
        if len(hz_idx_list) < 2:
            raise ValueError(f"Lysine {_resid} has <2 NZ hydrogens")

        nz_pos = protein_coords[nz_atom.molecule_atom_index]
        outward = nz_pos - centroid_vec
        outward_norm = np.linalg.norm(outward)
        if outward_norm < 1e-8:
            raise ValueError(f"Lysine {_resid} has zero outward vector")

        _site = ReactiveSite(
            resid=int(_resid),
            chain_id=str(nz_atom.metadata.get("chain_id", "A")),
            nz_index=nz_atom.molecule_atom_index,
            nz_position=nz_pos,
            ce_position=ce_pos,
            hz_indices=tuple(sorted(hz_idx_list)),
            outward_normal=outward / outward_norm,
        )
        reactive_sites[(_site.chain_id, _site.resid)] = _site
        logger.info(
            "  LYS %d: NZ=%d, %d Hz, outward=(%.2f, %.2f, %.2f)",
            _resid,
            _site.nz_index,
            len(hz_idx_list),
            *_site.outward_normal,
        )

    # ── Extract ff14SB charges for every protein atom ──────────────────────
    protein_ff = ForceField(PROTEIN_FF)
    protein_interchange = protein_ff.create_interchange(protein_top)
    _charges = protein_interchange["Electrostatics"].charges
    protein_charge_arr = np.array(
        [
            _charges[key].m_as(unit.elementary_charge)
            for key in sorted(_charges.keys(), key=lambda k: k.atom_indices[0])
        ]
    )
    logger.info(
        "Extracted ff14SB charges: %d atoms, total=%.4f e",
        len(protein_charge_arr),
        protein_charge_arr.sum(),
    )

    return (
        centroid_vec,
        protein_charge_arr,
        protein_coords,
        protein_mol,
        protein_top,
        reactive_sites,
    )


@app.cell
def _(
    AllChem,
    Chem,
    GENERATE_SCRIPT,
    MONOMER_GROUP_JSON,
    MONOMER_NAMES,
    Molecule,
    POLYMER_CACHE,
    POLYMER_SEQUENCE,
    POLYZYMD_PYTHON,
    PolymerPrototype,
    Quantity,
    RESIDUE_NAMES,
    json,
    logger,
    np,
    subprocess,
    unit,
):
    """Cell 3: Generate polymer via polyzymd subprocess, load PDB, identify NHS
    reactive group, compute product-state NAGL charges."""

    # ── Generate polymer PDB via subprocess ────────────────────────────────
    # The subprocess runs in the polyzymd build env (Python 3.11 + mbuild).
    _cmd = [
        POLYZYMD_PYTHON,
        GENERATE_SCRIPT,
        "--cache-dir", str(POLYMER_CACHE),
        "--monogrp-json", str(MONOMER_GROUP_JSON),
        "--sequence", POLYMER_SEQUENCE,
        "--monomer-names", json.dumps(MONOMER_NAMES),
        "--residue-names", json.dumps(RESIDUE_NAMES),
    ]

    # Check if PDB already cached
    _monomers_used = []
    for ch in POLYMER_SEQUENCE:
        if ch not in [m for m, _ in zip(MONOMER_NAMES.values(), range(len(_monomers_used)))]:
            _monomers_used.append(MONOMER_NAMES[ch])
    _unique_labels = []
    for ch in POLYMER_SEQUENCE:
        if ch not in _unique_labels:
            _unique_labels.append(ch)
    _prefix = "-".join(MONOMER_NAMES[ch] for ch in _unique_labels)
    _expected_pdb = POLYMER_CACHE / f"{_prefix}_seq={POLYMER_SEQUENCE}_{len(POLYMER_SEQUENCE)}-mer.pdb"

    if _expected_pdb.exists():
        logger.info("Using cached polymer PDB: %s", _expected_pdb)
        polymer_pdb_path = _expected_pdb
    else:
        logger.info("Generating polymer via subprocess...")
        _result = subprocess.run(
            _cmd, capture_output=True, text=True, timeout=300
        )
        if _result.returncode != 0:
            raise RuntimeError(
                f"Polymer generation failed:\nstdout={_result.stdout}\nstderr={_result.stderr}"
            )
        polymer_pdb_path = _result.stdout.strip().split("\n")[-1]
        logger.info("Generated polymer PDB: %s", polymer_pdb_path)

    # ── Load polymer PDB with RDKit ────────────────────────────────────────
    # RDKit reads CONECT records → full connectivity, but all bonds are single.
    # We keep this version for the final merge in Cell 6 (create_interchange
    # works fine with all-single-bond molecules).
    polymer_rdmol = Chem.MolFromPDBFile(str(polymer_pdb_path), removeHs=False, sanitize=True)
    if polymer_rdmol is None:
        raise ValueError(f"RDKit failed to load {polymer_pdb_path}")
    logger.info(
        "Polymer RDKit: %d atoms, %d bonds",
        polymer_rdmol.GetNumAtoms(),
        polymer_rdmol.GetNumBonds(),
    )

    # Create a SECOND copy with correct bond orders (C=O, S=O, C=C) for NAGL
    # charging.  DetermineBonds uses 3D geometry to infer bond orders.
    from rdkit.Chem import rdDetermineBonds
    _polymer_bo = Chem.MolFromPDBFile(str(polymer_pdb_path), removeHs=False, sanitize=False)
    rdDetermineBonds.DetermineConnectivity(_polymer_bo)
    rdDetermineBonds.DetermineBondOrders(_polymer_bo)
    Chem.SanitizeMol(_polymer_bo)
    logger.info("Applied DetermineBonds for NAGL charging model")

    # Also create an OpenFF Molecule (for conformer access and later merging)
    polymer_off = Molecule.from_rdkit(
        polymer_rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True,
    )

    # ── Identify NHS reactive carbon and leaving group ─────────────────────
    # In the PDB-loaded molecule all bonds are single, so SMARTS with "=O"
    # won't match.  Instead, use topology: the ester carbonyl C in the NH2
    # residue is the ONLY carbon bonded to exactly 2 oxygens.
    _nh2_indices = set()
    for _atom in polymer_rdmol.GetAtoms():
        _info = _atom.GetPDBResidueInfo()
        if _info and _info.GetResidueName().strip() == "NH2":
            _nh2_indices.add(_atom.GetIdx())

    if not _nh2_indices:
        raise ValueError("No NH2 (NHS) residue found in polymer PDB")

    _reactive_c = None
    _ester_o = None
    for _idx in _nh2_indices:
        _atom = polymer_rdmol.GetAtomWithIdx(_idx)
        if _atom.GetSymbol() != "C":
            continue
        _o_nbrs = [
            n.GetIdx() for n in _atom.GetNeighbors() if n.GetSymbol() == "O"
        ]
        if len(_o_nbrs) == 2:
            _reactive_c = _idx
            # The ester O is the one bonded to N (succinimide ring)
            for _o_idx in _o_nbrs:
                for _n in polymer_rdmol.GetAtomWithIdx(_o_idx).GetNeighbors():
                    if _n.GetSymbol() == "N" and _n.GetIdx() in _nh2_indices:
                        _ester_o = _o_idx
                        break
                if _ester_o is not None:
                    break
            break

    if _reactive_c is None or _ester_o is None:
        raise ValueError("Could not identify NHS ester carbonyl C and bridging O")

    # BFS from ester O (excluding reactive C) to find the leaving group
    _leaving: set[int] = set()
    _queue = [_ester_o]
    _visited = {_reactive_c}
    while _queue:
        _cur = _queue.pop(0)
        if _cur in _visited:
            continue
        _visited.add(_cur)
        _leaving.add(_cur)
        for _n in polymer_rdmol.GetAtomWithIdx(_cur).GetNeighbors():
            if _n.GetIdx() not in _visited:
                _queue.append(_n.GetIdx())

    _leaving_sorted = tuple(sorted(_leaving))
    _retained_sorted = tuple(
        i for i in range(polymer_rdmol.GetNumAtoms()) if i not in _leaving
    )
    logger.info(
        "NHS reactive C=%d, leaving group=%d atoms, retained=%d atoms",
        _reactive_c,
        len(_leaving_sorted),
        len(_retained_sorted),
    )

    # ── Compute product-state NAGL charges ─────────────────────────────────
    # Strategy: remove leaving group from the polymer, attach a small lysine-
    # sidechain cap (NH-CH2-CH2-CH2-CH2) at the reactive C to mimic the
    # product amide bond, then charge the whole thing with NAGL.
    # The cap charges are discarded; only retained-polymer charges are kept.

    # Step 1: Remove leaving group from the bond-order-corrected polymer mol
    # (using _polymer_bo so the product has correct C=O, S=O, C=C bond orders
    #  for accurate NAGL charge assignment).
    _core_rw = Chem.RWMol(_polymer_bo)
    for _idx in sorted(_leaving, reverse=True):
        _core_rw.RemoveAtom(_idx)
    _core = _core_rw.GetMol()
    Chem.SanitizeMol(_core)

    _old_to_core: dict[int, int] = {}
    _cursor = 0
    for _old_idx in range(polymer_rdmol.GetNumAtoms()):
        if _old_idx not in _leaving:
            _old_to_core[_old_idx] = _cursor
            _cursor += 1
    _rc_core = _old_to_core[_reactive_c]

    # Step 2: Build lysine sidechain cap: NH-CH2-CH2-CH2-CH2
    # from_smiles("NCCCC") gives NH2-CH2-CH2-CH2-CH3 (primary amine, 3 bonds on N).
    # We must remove one H from N before bonding to reactive C (amide = 3 bonds on N).
    _cap = Molecule.from_smiles("NCCCC", allow_undefined_stereo=True)
    _cap.generate_conformers(n_conformers=1)
    _cap_rd = _cap.to_rdkit()
    _cap_rw = Chem.RWMol(_cap_rd)
    _cap_n = _cap_rw.GetAtomWithIdx(0)
    _h_drop = None
    for _b in _cap_n.GetBonds():
        _nbr = _b.GetOtherAtomIdx(0)
        if _cap_rw.GetAtomWithIdx(_nbr).GetSymbol() == "H":
            _h_drop = _nbr
            break
    if _h_drop is not None:
        _cap_rw.RemoveAtom(_h_drop)
    _cap_rd = _cap_rw.GetMol()
    Chem.SanitizeMol(_cap_rd)

    # Step 3: Combine core + cap, add amide bond
    _merged_for_charge = AllChem.CombineMols(_core, _cap_rd)
    _merged_rw = Chem.RWMol(_merged_for_charge)
    _cap_offset = _core.GetNumAtoms()
    _merged_rw.AddBond(_rc_core, _cap_offset, Chem.BondType.SINGLE)
    _product_rd = _merged_rw.GetMol()
    Chem.SanitizeMol(_product_rd)

    _product_off = Molecule.from_rdkit(
        _product_rd, allow_undefined_stereo=True, hydrogens_are_explicit=True,
    )
    # No conformer needed — NAGL charges from graph, not geometry.

    # Step 4: Assign NAGL charges
    try:
        from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        _product_off.assign_partial_charges(
            partial_charge_method="openff-gnn-am1bcc-0.1.0-rc.3.pt",
            toolkit_registry=NAGLToolkitWrapper(),
        )
        logger.info("Product-state charges assigned via NAGL")
    except (ImportError, ValueError, RuntimeError) as _exc:
        logger.warning("NAGL failed (%s), falling back to am1bcc", _exc)
        _product_off.assign_partial_charges(partial_charge_method="am1bcc")

    # Step 5: Extract charges for retained polymer atoms only
    _product_charges_all = np.array(
        [c.m_as(unit.elementary_charge) for c in _product_off.partial_charges]
    )
    # The first _core.GetNumAtoms() entries in _product_off correspond to core atoms
    # which map 1:1 with retained polymer atoms.
    _product_charges_core = _product_charges_all[: _core.GetNumAtoms()]

    # Map core indices back to original polymer indices for the retained atoms
    _core_to_old = {v: k for k, v in _old_to_core.items()}
    _retained_charges = np.array(
        [_product_charges_core[_old_to_core[i]] for i in _retained_sorted]
    )
    logger.info(
        "Product charges: %d retained atoms, sum=%.6f e",
        len(_retained_charges),
        _retained_charges.sum(),
    )

    # Also extract cap-atom charges for patching the conjugated lysine.
    # Cap atoms in the product molecule start at index _cap_offset.
    # Cap structure after H removal: N(0) - C(1) - C(2) - C(3) - C(4)
    # These map to lysine atoms:     NZ     CE     CD     CG     CB
    _cap_charge_map = {}
    _cap_atom_names = ["NZ", "CE", "CD", "CG", "CB"]
    _cap_heavy_cursor = 0
    for _ci in range(_cap_offset, _product_off.n_atoms):
        _ca = _product_off.atom(_ci)
        if _ca.atomic_number > 1:
            if _cap_heavy_cursor < len(_cap_atom_names):
                _cap_charge_map[_cap_atom_names[_cap_heavy_cursor]] = float(
                    _product_charges_all[_ci]
                )
            _cap_heavy_cursor += 1
        else:
            # Cap H atoms — the first one bonded to N is HZ1
            if _cap_heavy_cursor == 1:  # just finished N, this H is on N
                _cap_charge_map["HZ1"] = float(_product_charges_all[_ci])

    logger.info("Cap charge map (for LYS patch): %s", _cap_charge_map)

    polymer_prototype = PolymerPrototype(
        off_molecule=polymer_off,
        rdmol=polymer_rdmol,
        reactive_carbon_idx=_reactive_c,
        leaving_group_indices=_leaving_sorted,
        retained_atom_indices=_retained_sorted,
        product_charges=_retained_charges,
    )

    cap_charge_map = _cap_charge_map
    return cap_charge_map, polymer_prototype


@app.cell
def _(
    N_CONJUGATED,
    RNG_SEED,
    logger,
    np,
    polymer_prototype,
    reactive_sites,
):
    """Cell 4: Assign polymer prototypes to reactive lysine sites."""

    site_list = sorted(reactive_sites.values(), key=lambda s: s.resid)
    if len(site_list) < N_CONJUGATED:
        raise ValueError(
            f"Requested {N_CONJUGATED} conjugations but only {len(site_list)} sites"
        )

    rng = np.random.default_rng(RNG_SEED)
    site_indices = rng.choice(len(site_list), size=N_CONJUGATED, replace=False)

    # All conjugation sites use the same polymer prototype in this PoC.
    assignments: list[tuple] = [
        (site_list[int(idx)], polymer_prototype) for idx in site_indices
    ]
    for _site, _poly in assignments:
        logger.info("  LYS %d <- polymer (%d atoms)", _site.resid, _poly.off_molecule.n_atoms)
    return (assignments,)


@app.cell
def _(
    BOND_ANGLE_TOLERANCE,
    COLLISION_DISTANCE,
    KDTree,
    N_CONFORMERS,
    N_CONE_SAMPLES,
    N_SPIN_SAMPLES,
    PlacementResult,
    RNG_SEED,
    TARGET_BOND_ANGLE,
    TARGET_BOND_LENGTH,
    assignments,
    centroid_vec,
    logger,
    np,
    protein_coords,
    protein_mol,
    reactive_sites,
):
    """Cell 5: Geometric placement of polymers at lysine sites.

    Uses cone + spin sampling around the CE→NZ direction to place each
    polymer's reactive carbon at the target amide bond length from NZ,
    then spins the polymer around the NZ→C axis to find a collision-free
    orientation.
    """

    def _rodrigues_rotate(v: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
        ct = np.cos(theta)
        st = np.sin(theta)
        return v * ct + np.cross(axis, v) * st + axis * np.dot(axis, v) * (1 - ct)

    def _place_one(
        site,
        polymer,
        occupied_coords: np.ndarray,
        centroid: np.ndarray,
    ) -> PlacementResult:
        rng_local = np.random.default_rng(RNG_SEED + int(site.resid))

        ce_to_nz = site.nz_position - site.ce_position
        ce_to_nz_hat = ce_to_nz / np.linalg.norm(ce_to_nz)

        nz_outward = site.nz_position - centroid
        proj = np.dot(nz_outward, ce_to_nz_hat) * ce_to_nz_hat
        perp = nz_outward - proj
        perp_norm = np.linalg.norm(perp)
        if perp_norm < 1e-8:
            arb = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(arb, ce_to_nz_hat)) > 0.9:
                arb = np.array([0.0, 1.0, 0.0])
            perp = arb - np.dot(arb, ce_to_nz_hat) * ce_to_nz_hat
            perp_norm = np.linalg.norm(perp)
        outward_perp = perp / perp_norm

        cone_half_angle = np.radians(180.0 - TARGET_BOND_ANGLE)

        retained_idx = np.asarray(polymer.retained_atom_indices, dtype=int)
        collision_mask = retained_idx != polymer.reactive_carbon_idx
        collision_checked_idx = retained_idx[collision_mask]

        occupied_tree = KDTree(occupied_coords)

        # Generate conformers from the precursor molecule
        polymer_mol_copy = Molecule(polymer.off_molecule)
        try:
            polymer_mol_copy.generate_conformers(n_conformers=N_CONFORMERS)
        except Exception:
            pass
        conformers = [c.m_as("angstrom") for c in polymer_mol_copy.conformers]

        best_coords = None
        best_clearance = -1.0
        total_attempts = 0

        for conf_coords in conformers:
            reactive_pos = conf_coords[polymer.reactive_carbon_idx]
            retained_com = conf_coords[retained_idx].mean(axis=0)
            com_dir = retained_com - reactive_pos
            com_norm = np.linalg.norm(com_dir)
            if com_norm < 1e-8:
                continue
            com_hat = com_dir / com_norm

            for _cone_idx in range(N_CONE_SAMPLES):
                spin_axis_angle = rng_local.uniform(0, 2 * np.pi)
                spun_perp = _rodrigues_rotate(outward_perp, ce_to_nz_hat, spin_axis_angle)
                nz_to_c_dir = _rodrigues_rotate(ce_to_nz_hat, spun_perp, cone_half_angle)
                nz_to_c_dir = nz_to_c_dir / np.linalg.norm(nz_to_c_dir)

                c_target = site.nz_position + TARGET_BOND_LENGTH * nz_to_c_dir

                cross = np.cross(com_hat, nz_to_c_dir)
                cross_norm = np.linalg.norm(cross)
                dot = np.dot(com_hat, nz_to_c_dir)
                if cross_norm < 1e-8:
                    if dot > 0:
                        R_align = np.eye(3)
                    else:
                        arb2 = np.array([1.0, 0.0, 0.0])
                        if abs(np.dot(arb2, com_hat)) > 0.9:
                            arb2 = np.array([0.0, 1.0, 0.0])
                        flip_axis = np.cross(com_hat, arb2)
                        flip_axis /= np.linalg.norm(flip_axis)
                        R_align = 2 * np.outer(flip_axis, flip_axis) - np.eye(3)
                else:
                    axis = cross / cross_norm
                    angle = np.arctan2(cross_norm, dot)
                    K = np.array([
                        [0, -axis[2], axis[1]],
                        [axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0],
                    ])
                    R_align = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

                centered = conf_coords - reactive_pos
                aligned = (R_align @ centered.T).T + c_target

                for _spin_idx in range(N_SPIN_SAMPLES):
                    total_attempts += 1
                    spin_angle = rng_local.uniform(0, 2 * np.pi)
                    ct = np.cos(spin_angle)
                    st = np.sin(spin_angle)
                    ux, uy, uz = nz_to_c_dir
                    R_spin = np.array([
                        [ct + ux*ux*(1-ct), ux*uy*(1-ct) - uz*st, ux*uz*(1-ct) + uy*st],
                        [uy*ux*(1-ct) + uz*st, ct + uy*uy*(1-ct), uy*uz*(1-ct) - ux*st],
                        [uz*ux*(1-ct) - uy*st, uz*uy*(1-ct) + ux*st, ct + uz*uz*(1-ct)],
                    ])
                    spun = (R_spin @ (aligned - c_target).T).T + c_target

                    check_pts = spun[collision_checked_idx]
                    dists, _ = occupied_tree.query(check_pts)
                    min_dist = float(np.min(dists))
                    if min_dist < COLLISION_DISTANCE:
                        continue

                    nz_to_ce = site.ce_position - site.nz_position
                    nz_to_c = spun[polymer.reactive_carbon_idx] - site.nz_position
                    denom = np.linalg.norm(nz_to_ce) * np.linalg.norm(nz_to_c)
                    if denom == 0.0:
                        continue
                    cosine = np.dot(nz_to_ce, nz_to_c) / denom
                    bond_angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
                    if abs(bond_angle - TARGET_BOND_ANGLE) > BOND_ANGLE_TOLERANCE:
                        continue

                    if min_dist > best_clearance:
                        best_coords = spun.copy()
                        best_clearance = min_dist

        if best_coords is None:
            raise RuntimeError(
                f"No valid placement for LYS {site.resid} after {total_attempts} attempts"
            )
        logger.info(
            "LYS %d: clearance %.2f A (%d attempts)",
            site.resid, best_clearance, total_attempts,
        )
        return PlacementResult(
            site=site,
            polymer=polymer,
            placed_coords=best_coords,
            min_protein_distance=best_clearance,
        )

    # ── Exclude only the HZ atoms being removed + NZ from occupied coords ──
    # NZ is excluded because the reactive C will bond directly to it at close
    # range.  All other lysine atoms (CE, CD, CG, CB, backbone) remain as
    # collision targets so the polymer doesn't clash with them.
    exclude_indices: set[int] = set()
    for _site_key in reactive_sites:
        _site_obj = reactive_sites[_site_key]
        # Exclude NZ
        exclude_indices.add(_site_obj.nz_index)
        # Exclude the HZ that will be removed (last 2 of the sorted list)
        for _hz in _site_obj.hz_indices[1:]:
            exclude_indices.add(_hz)

    keep_mask = np.ones(len(protein_coords), dtype=bool)
    keep_mask[list(exclude_indices)] = False
    occupied = protein_coords[keep_mask].copy()

    placement_results: list[PlacementResult] = []
    for _site, _polymer in assignments:
        _placement = _place_one(_site, _polymer, occupied, centroid_vec)
        placement_results.append(_placement)
        retained_coords = _placement.placed_coords[
            np.asarray(_polymer.retained_atom_indices, dtype=int)
        ]
        occupied = np.vstack((occupied, retained_coords))
        logger.info(
            "Placed polymer at LYS %d with clearance %.3f A",
            _site.resid, _placement.min_protein_distance,
        )

    return (placement_results,)


@app.cell
def _(
    Chem,
    Molecule,
    Quantity,
    cap_charge_map,
    logger,
    np,
    placement_results,
    protein_charge_arr,
    protein_coords,
    protein_mol,
    unit,
):
    """Cell 6: RDKit graph surgery — merge protein + polymers into a single
    molecule with crosslink bonds.  Build the combined charge vector.

    Strategy
    --------
    1. Start with the protein RDKit mol.
    2. For each conjugation site:
       a. Remove 2 of 3 HZ from the lysine NZ (keep one for amide N-H).
       b. Remove the NHS leaving group from the polymer.
       c. CombineMols the growing molecule with the trimmed polymer.
       d. AddBond between NZ and the reactive carbonyl C.
    3. Sanitize the merged molecule.
    4. Build the coordinate array by combining protein + placed-polymer coords.
    5. Build the charge vector: ff14SB for protein, NAGL product-state for
       polymer, with cap-derived patches for the modified lysine atoms near
       the junction.
    """

    # Start from protein RDKit mol
    _prot_rd = Chem.MolFromPDBFile(
        protein_mol.conformers[0]._parent.to_topology().to_pdb()
        if False  # placeholder — we need the actual PDB path
        else "/home/joelaforet/Shirts-Lab-Linux/openff-pablo/SBMA_EGMA_50_50_water_363K/structures/NH3_terminal_His_proton_updated.pdb",
        removeHs=False,
        sanitize=True,
    )
    if _prot_rd is None:
        raise ValueError("Failed to load protein with RDKit")

    # We need to map between the OpenFF protein_mol atom ordering and
    # the RDKit mol ordering.  Topology.from_pdb and RDKit MolFromPDBFile
    # should give the same atom order when both read the same PDB.

    # ── Track cumulative modifications ─────────────────────────────────────
    # As we remove H atoms and add polymers, indices shift.  Track via
    # old→new maps at each stage.
    _current_rd = _prot_rd
    _n_prot_orig = _prot_rd.GetNumAtoms()

    # Global map: original-protein-idx → current-merged-idx
    _prot_global_map: dict[int, int] = {i: i for i in range(_n_prot_orig)}
    _prot_removed: set[int] = set()  # original protein indices removed so far

    # For charge building: track which atoms are protein vs polymer
    # and store per-atom charges
    _charges: list[float] = []
    _source_labels: list[str] = []  # "protein" or "polymer_N"

    # For coordinate building
    _coord_blocks: list[np.ndarray] = [protein_coords.copy()]
    _prot_h_removed_indices: list[int] = []  # original protein indices of removed H

    # Collect per-site info for PDB assembly
    _polymer_atom_ranges: list[tuple[int, int]] = []  # (start, end) in merged mol

    for _pi, _placement in enumerate(placement_results):
        _site = _placement.site
        _polymer = _placement.polymer

        # ── Identify which 2 HZ to remove from this lysine ────────────────
        # Keep the first HZ, remove the other two.
        # hz_indices are in the ORIGINAL protein_mol ordering.
        _hz_to_remove = list(_site.hz_indices[1:])  # remove last 2, keep first
        _prot_removed.update(_hz_to_remove)
        _prot_h_removed_indices.extend(_hz_to_remove)

    # ── Build modified protein (remove all HZ at once) ─────────────────────
    _prot_rw = Chem.RWMol(_prot_rd)
    for _idx in sorted(_prot_removed, reverse=True):
        _prot_rw.RemoveAtom(_idx)
    _prot_mod = _prot_rw.GetMol()
    Chem.SanitizeMol(_prot_mod)

    # Build old→new map for protein atoms after H removal
    _prot_old_to_new: dict[int, int] = {}
    _cursor = 0
    for _old in range(_n_prot_orig):
        if _old not in _prot_removed:
            _prot_old_to_new[_old] = _cursor
            _cursor += 1
    _n_prot_mod = _prot_mod.GetNumAtoms()

    # Build protein charge array (remove HZ entries, redistribute charge)
    _prot_charges_mod: list[float] = []
    _prot_keep_mask = np.ones(_n_prot_orig, dtype=bool)
    _prot_keep_mask[list(_prot_removed)] = False
    for _old in range(_n_prot_orig):
        if _old not in _prot_removed:
            _prot_charges_mod.append(float(protein_charge_arr[_old]))

    # Build protein coordinate array (remove H entries)
    _prot_coords_mod = protein_coords[_prot_keep_mask]

    # ── Now merge polymers one by one ──────────────────────────────────────
    _merged = _prot_mod
    _merged_charges = list(_prot_charges_mod)
    _merged_coords = _prot_coords_mod.copy()

    for _pi, _placement in enumerate(placement_results):
        _site = _placement.site
        _polymer = _placement.polymer
        _nz_new = _prot_old_to_new[_site.nz_index]

        # Remove leaving group from polymer
        _poly_rd = _polymer.rdmol
        _poly_rw = Chem.RWMol(_poly_rd)
        for _idx in sorted(_polymer.leaving_group_indices, reverse=True):
            _poly_rw.RemoveAtom(_idx)
        _poly_mod = _poly_rw.GetMol()
        Chem.SanitizeMol(_poly_mod)

        _poly_old_to_new: dict[int, int] = {}
        _cur = 0
        for _old in range(_poly_rd.GetNumAtoms()):
            if _old not in set(_polymer.leaving_group_indices):
                _poly_old_to_new[_old] = _cur
                _cur += 1
        _rc_new = _poly_old_to_new[_polymer.reactive_carbon_idx]

        # Merge
        _offset = _merged.GetNumAtoms()
        _combined = Chem.CombineMols(_merged, _poly_mod)
        _combined_rw = Chem.RWMol(_combined)
        _combined_rw.AddBond(_nz_new, _offset + _rc_new, Chem.BondType.SINGLE)
        _merged = _combined_rw.GetMol()
        Chem.SanitizeMol(_merged)

        # Append polymer charges (retained atoms only, in retained order)
        for _ri, _orig_idx in enumerate(_polymer.retained_atom_indices):
            _merged_charges.append(float(_polymer.product_charges[_ri]))

        # Append polymer coordinates (retained atoms, placed)
        _poly_placed = _placement.placed_coords
        _retained_coords = _poly_placed[np.array(_polymer.retained_atom_indices)]
        _merged_coords = np.vstack([_merged_coords, _retained_coords])

        logger.info(
            "Merged polymer %d at LYS %d: merged total=%d atoms",
            _pi, _site.resid, _merged.GetNumAtoms(),
        )

    # ── Patch lysine junction charges ──────────────────────────────────────
    # For each conjugated lysine, replace the ff14SB charges on
    # NZ, CE, CD, CG, CB, HZ1 with the cap-derived charges.
    for _placement in placement_results:
        _site = _placement.site
        for _atom_name, _cap_q in cap_charge_map.items():
            # Find the protein atom by name in the original ordering
            for _atom in protein_mol.atoms:
                if (
                    _atom.metadata.get("residue_number") == str(_site.resid)
                    and _atom.metadata.get("chain_id") == _site.chain_id
                    and _atom.metadata.get("residue_name") == "LYS"
                    and _atom.name == _atom_name
                ):
                    _orig_idx = _atom.molecule_atom_index
                    if _orig_idx in _prot_old_to_new:
                        _new_idx = _prot_old_to_new[_orig_idx]
                        _old_q = _merged_charges[_new_idx]
                        _merged_charges[_new_idx] = _cap_q
                        logger.info(
                            "  Patched LYS %d %s: %.4f → %.4f",
                            _site.resid, _atom_name, _old_q, _cap_q,
                        )
                    break

    # ── Convert to OpenFF Molecule ─────────────────────────────────────────
    # Remove RDKit conformers first — they carry stale coordinates from the
    # original PDB loads (pre-placement).  We set the correct coordinates
    # from _merged_coords below.
    _merged_noconf = Chem.RWMol(_merged)
    _merged_noconf.RemoveAllConformers()
    _merged_noconf = _merged_noconf.GetMol()

    conjugate_off = Molecule.from_rdkit(
        _merged_noconf, allow_undefined_stereo=True, hydrogens_are_explicit=True,
    )
    conjugate_off.add_conformer(Quantity(_merged_coords, unit.angstrom))

    # ── Fix formal charges on protein atoms ────────────────────────────────
    # RDKit assigns different formal charges than OpenFF/CCD for some residues
    # (e.g. ARG NH1 gets FC=+1 in RDKit but FC=0 in OpenFF).  We correct
    # these where safe (i.e., where the new FC doesn't violate valence rules).
    #
    # Conjugated lysine NZ changes from +1 (protonated NH3+) to
    # 0 (amide N) after removing 2 HZ and forming the crosslink bond.
    #
    # We do NOT change HIS ND1 because it is aromatic and FC=0 would break
    # kekulization.
    _conjugated_nz_indices = {
        _prot_old_to_new[p.site.nz_index] for p in placement_results
    }
    _n_fc_fixed = 0
    for _orig_idx in range(protein_mol.n_atoms):
        if _orig_idx in _prot_removed:
            continue
        _new_idx = _prot_old_to_new[_orig_idx]
        if _new_idx in _conjugated_nz_indices:
            _target_fc = 0
        else:
            _target_fc = protein_mol.atom(_orig_idx).formal_charge.m
        _current_fc = conjugate_off.atom(_new_idx).formal_charge.m
        if _target_fc == _current_fc:
            continue
        # Safety check: don't change FC on aromatic atoms (e.g. HIS ND1)
        # which could break kekulization
        _rd_atom = _merged.GetAtomWithIdx(_new_idx)
        if _rd_atom.GetIsAromatic():
            logger.info(
                "  Skipping FC fix for aromatic atom %d (%s): keeping FC=%d",
                _orig_idx,
                protein_mol.atom(_orig_idx).name,
                _current_fc,
            )
            continue
        conjugate_off.atom(_new_idx).formal_charge = _target_fc
        _n_fc_fixed += 1
        logger.info(
            "  Fixed formal charge: protein atom %d (%s) %d → %d",
            _orig_idx,
            protein_mol.atom(_orig_idx).name,
            _current_fc,
            _target_fc,
        )
    logger.info("Fixed %d formal charges", _n_fc_fixed)

    # Assign charges
    _charge_arr = np.array(_merged_charges)
    # Ensure total charge matches formal charge
    _formal_charge = sum(a.formal_charge.m for a in conjugate_off.atoms)
    _charge_sum = _charge_arr.sum()
    _charge_diff = _formal_charge - _charge_sum
    if abs(_charge_diff) > 0.01:
        logger.warning(
            "Charge mismatch: formal=%.4f, sum=%.4f, diff=%.4f — redistributing",
            _formal_charge, _charge_sum, _charge_diff,
        )
        # Spread the difference evenly across all atoms
        _charge_arr += _charge_diff / len(_charge_arr)
    conjugate_off.partial_charges = Quantity(_charge_arr, unit.elementary_charge)

    logger.info(
        "Conjugate molecule: %d atoms, formal charge=%.0f, partial charge sum=%.4f",
        conjugate_off.n_atoms,
        _formal_charge,
        _charge_arr.sum(),
    )

    merged_rdmol = _merged
    prot_old_to_new = _prot_old_to_new
    prot_removed_indices = _prot_removed
    return (conjugate_off, merged_rdmol, prot_old_to_new, prot_removed_indices)


@app.cell
def _(
    ASSEMBLED_PDB,
    Path,
    logger,
    np,
    placement_results,
    protein_coords,
    protein_mol,
    prot_removed_indices,
    write_conect_line,
    write_pdb_line,
):
    """Cell 7: Write assembled PDB for visualization and minimization."""

    # Parse original protein PDB for ATOM lines
    _pdb_path = (
        "/home/joelaforet/Shirts-Lab-Linux/openff-pablo/SBMA_EGMA_50_50_water_363K/"
        "structures/NH3_terminal_His_proton_updated.pdb"
    )
    _atom_records: list[dict] = []
    with open(_pdb_path, "r", encoding="utf-8") as _f:
        for _line in _f:
            if _line.startswith(("ATOM", "HETATM")):
                _atom_records.append({
                    "line": _line,
                    "serial": int(_line[6:11]),
                    "name": _line[12:16].strip(),
                    "resname": _line[17:20].strip(),
                    "chain": _line[21].strip() or "A",
                    "resid": int(_line[22:26]),
                })

    # Map original protein atom index to PDB serial
    # (assumes PDB ATOM lines are in the same order as Topology.from_pdb atoms)
    _orig_idx_to_serial: dict[int, int] = {}
    for _i, _rec in enumerate(_atom_records):
        _orig_idx_to_serial[_i] = _rec["serial"]

    # Write protein atoms, skipping removed H
    kept_lines: list[str] = []
    _serial_lookup: dict[int, int] = {}  # original idx → serial in output
    _existing_serials: set[int] = set()
    for _orig_idx, _rec in enumerate(_atom_records):
        if _orig_idx in prot_removed_indices:
            continue
        kept_lines.append(_rec["line"])
        _serial_lookup[_orig_idx] = _rec["serial"]
        _existing_serials.add(_rec["serial"])

    _next_serial = max(_existing_serials) + 1

    # Write polymer atoms
    polymer_lines: list[str] = []
    conect_map: dict[int, set[int]] = {}

    for _pi, _placement in enumerate(placement_results):
        _polymer = _placement.polymer
        _poly_chain = "C"
        _poly_resid = 1001 + _pi

        _serial_map: dict[int, int] = {}  # polymer retained idx → serial
        for _ri, _orig_idx in enumerate(_polymer.retained_atom_indices):
            _atom = _polymer.off_molecule.atoms[_orig_idx]
            _x, _y, _z = _placement.placed_coords[_orig_idx]
            # Use PDB residue info from RDKit if available
            _rd_atom = _polymer.rdmol.GetAtomWithIdx(_orig_idx)
            _rd_info = _rd_atom.GetPDBResidueInfo()
            _atom_name = _rd_info.GetName().strip() if _rd_info else f"{_atom.symbol}{_ri:03d}"
            _res_name = _rd_info.GetResidueName().strip() if _rd_info else "POL"

            polymer_lines.append(
                write_pdb_line(
                    serial=_next_serial,
                    atom_name=_atom_name,
                    residue_name=_res_name,
                    chain_id=_poly_chain,
                    residue_number=_rd_info.GetResidueNumber() + 1000 * (_pi + 1) if _rd_info else _poly_resid,
                    x=float(_x), y=float(_y), z=float(_z),
                    element=_atom.symbol,
                    record="HETATM",
                )
            )
            _serial_map[_orig_idx] = _next_serial
            _next_serial += 1

        # CONECT records for intra-polymer bonds (retained atoms only)
        _retained_set = set(_polymer.retained_atom_indices)
        for _bond in _polymer.off_molecule.bonds:
            _i = _bond.atom1_index
            _j = _bond.atom2_index
            if _i in _retained_set and _j in _retained_set:
                _si = _serial_map[_i]
                _sj = _serial_map[_j]
                conect_map.setdefault(_si, set()).add(_sj)
                conect_map.setdefault(_sj, set()).add(_si)

        # CONECT for the crosslink bond (NZ — reactive C)
        _nz_serial = _serial_lookup[_placement.site.nz_index]
        _rc_serial = _serial_map[_polymer.reactive_carbon_idx]
        conect_map.setdefault(_nz_serial, set()).add(_rc_serial)
        conect_map.setdefault(_rc_serial, set()).add(_nz_serial)

    # Write the assembled PDB
    assembled_path = Path(ASSEMBLED_PDB)
    with open(assembled_path, "w", encoding="utf-8") as _f:
        for _line in kept_lines:
            _f.write(_line)
        for _line in polymer_lines:
            _f.write(_line)
        _f.write("TER\n")
        for _serial in sorted(conect_map):
            _bonded = sorted(conect_map[_serial])
            _f.write(write_conect_line(_serial, _bonded))
        _f.write("END\n")

    logger.info("Wrote assembled PDB: %s", assembled_path)
    return (assembled_path,)


@app.cell
def _(
    ForceField,
    PROTEIN_FF,
    SMALL_MOL_FF,
    conjugate_off,
    logger,
):
    """Cell 8: Create Interchange from the merged conjugate molecule."""

    conjugate_top = conjugate_off.to_topology()

    ff = ForceField(PROTEIN_FF, SMALL_MOL_FF)
    interchange = ff.create_interchange(
        conjugate_top, charge_from_molecules=[conjugate_off]
    )
    logger.info(
        "Created interchange: %d atoms, %d bonds",
        interchange.topology.n_atoms,
        interchange.topology.n_bonds,
    )
    return (interchange,)


@app.cell
def _(
    LINKAGE_NEIGHBORHOOD_BONDS,
    MINIMIZATION_MAX_ITER,
    MINIMIZED_PDB,
    Molecule,
    assembled_path,
    conjugate_off,
    interchange,
    logger,
    placement_results,
    unit,
):
    """Cell 9: OpenMM local energy minimization with position restraints.

    Only atoms within LINKAGE_NEIGHBORHOOD_BONDS graph bonds of the
    crosslink (NZ and reactive C) are allowed to move.  All other atoms
    are restrained to their initial positions.
    """
    try:
        import openmm
        from openmm import app as openmm_app
        from openmm import unit as omm_unit
    except ImportError as _exc:
        raise ImportError("OpenMM required for minimization") from _exc

    def _atoms_within_n_bonds(molecule: Molecule, seeds: set[int], n: int) -> set[int]:
        reached = set(seeds)
        frontier = set(seeds)
        for _ in range(n):
            next_frontier: set[int] = set()
            for idx in frontier:
                for bond in molecule.atom(idx).bonds:
                    nbr = (
                        bond.atom1_index if bond.atom2_index == idx else bond.atom2_index
                    )
                    if nbr not in reached:
                        reached.add(nbr)
                        next_frontier.add(nbr)
            frontier = next_frontier
        return reached

    # Find linkage seed atoms (NZ and reactive C in the merged molecule).
    # In the merged molecule, protein atoms come first (minus removed H),
    # then polymer retained atoms appended per placement.
    _n_prot_atoms = conjugate_off.n_atoms - sum(
        len(p.polymer.retained_atom_indices) for p in placement_results
    )

    linkage_seeds: set[int] = set()
    # For each placement, the NZ index in the merged mol is the same as
    # in _prot_old_to_new (passed via conjugate_off construction).
    # The reactive C is at offset + _poly_old_to_new[reactive_carbon_idx].
    # Rather than tracking all those maps, we find them by looking at
    # the bond between a protein-range atom and a polymer-range atom.
    for _bond in conjugate_off.bonds:
        _i = _bond.atom1_index
        _j = _bond.atom2_index
        if (_i < _n_prot_atoms) != (_j < _n_prot_atoms):
            # Cross-boundary bond = crosslink
            linkage_seeds.add(_i)
            linkage_seeds.add(_j)

    if not linkage_seeds:
        raise ValueError("No crosslink bonds found in conjugate molecule")
    logger.info("Linkage seeds: %s", linkage_seeds)

    mobile = _atoms_within_n_bonds(
        conjugate_off, linkage_seeds, LINKAGE_NEIGHBORHOOD_BONDS
    )
    logger.info("Mobile atoms (within %d bonds): %d", LINKAGE_NEIGHBORHOOD_BONDS, len(mobile))

    system = interchange.to_openmm_system()
    topology_omm = interchange.to_openmm_topology()
    _coords = conjugate_off.conformers[0].m_as(unit.angstrom)

    restraint = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter(
        "k", 1000.0 * omm_unit.kilojoule_per_mole / omm_unit.nanometer**2
    )
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    coords_nm = _coords * 0.1
    for _idx in range(conjugate_off.n_atoms):
        if _idx in mobile:
            continue
        x0, y0, z0 = coords_nm[_idx]
        restraint.addParticle(_idx, [float(x0), float(y0), float(z0)])
    system.addForce(restraint)

    integrator = openmm.VerletIntegrator(0.001 * omm_unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("Reference")
    simulation = openmm_app.Simulation(topology_omm, system, integrator, platform)
    simulation.context.setPositions(coords_nm * omm_unit.nanometers)

    before = simulation.context.getState(getEnergy=True)
    energy_before = before.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    openmm.LocalEnergyMinimizer.minimize(
        simulation.context, tolerance=10.0, maxIterations=MINIMIZATION_MAX_ITER
    )
    after = simulation.context.getState(getEnergy=True, getPositions=True)
    energy_after = after.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
    minimized_coords = after.getPositions(asNumpy=True).value_in_unit(omm_unit.angstrom)

    # Update conformer
    for _idx, xyz in enumerate(minimized_coords):
        conjugate_off.conformers[0][_idx, :] = xyz * unit.angstrom

    # Write minimized PDB by updating coordinates in the assembled PDB
    with open(assembled_path, "r", encoding="utf-8") as _f:
        _old_lines = _f.readlines()
    _atom_lines = [l for l in _old_lines if l.startswith(("ATOM", "HETATM"))]
    _other_lines = [l for l in _old_lines if not l.startswith(("ATOM", "HETATM"))]

    _updated: list[str] = []
    for _i, _line in enumerate(_atom_lines):
        _x, _y, _z = minimized_coords[_i]
        _updated.append(f"{_line[:30]}{_x:8.3f}{_y:8.3f}{_z:8.3f}{_line[54:]}")

    output_path = MINIMIZED_PDB
    with open(output_path, "w", encoding="utf-8") as _f:
        for _line in _updated:
            _f.write(_line)
        for _line in _other_lines:
            _f.write(_line)

    logger.info(
        "Minimization: E_before=%.2f, E_after=%.2f kJ/mol",
        energy_before, energy_after,
    )
    return energy_after, energy_before, output_path


@app.cell
def _(
    energy_after,
    energy_before,
    logger,
    output_path,
    placement_results,
):
    """Cell 10: Summary and visualization."""
    import marimo as mo

    summary_rows = []
    for _i, _p in enumerate(placement_results):
        summary_rows.append({
            "id": _i,
            "resid": _p.site.resid,
            "clearance_A": round(_p.min_protein_distance, 3),
        })

    logger.info("Summary: %s", summary_rows)

    try:
        import py3Dmol

        with open(output_path, "r", encoding="utf-8") as _f:
            pdb_block = _f.read()
        viewer = py3Dmol.view(width=900, height=600)
        viewer.addModel(pdb_block, "pdb")
        viewer.setStyle({"chain": "A"}, {"cartoon": {"color": "spectrum"}})
        _colors = ["greenCarbon", "cyanCarbon", "orangeCarbon", "magentaCarbon"]
        for _i, _p in enumerate(placement_results):
            viewer.addStyle(
                {"chain": "C"},
                {"stick": {"colorscheme": _colors[_i % len(_colors)]}},
            )
        viewer.zoomTo()
        viz = mo.Html(viewer._make_html())
    except ImportError:
        viz = mo.md(f"py3Dmol not available. PDB: `{output_path}`")

    summary = {
        "n_conjugates": len(placement_results),
        "energy_before_kj_mol": float(energy_before),
        "energy_after_kj_mol": float(energy_after),
        "output_pdb": str(output_path),
        "rows": summary_rows,
    }
    return summary, viz


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
