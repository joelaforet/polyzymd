import marimo

__generated_with = "0.22.5"
app = marimo.App()


@app.cell
def _():
    """Cell 1: Configuration, imports, and shared utilities."""
    import logging
    from dataclasses import dataclass
    from os import environ
    from pathlib import Path

    import numpy as np
    from openff.toolkit import ForceField, Molecule, Topology
    from openff.units import Quantity, unit
    from polymerist.polymers.monomers import MonomerGroup
    from rdkit import Chem
    from scipy.spatial import KDTree
    from scipy.spatial.transform import Rotation

    from polyzymd.builders.fragment_generator import FragmentGenerator
    from polyzymd.builders.polymer_generator import PolymerGenerator

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("conjugation_poc")

    # ── Paths ──────────────────────────────────────────────────────────────
    _POC_DIR = Path(__file__).resolve().parent

    PROTEIN_PDB = environ.get(
        "CONJUGATION_PROTEIN_PDB",
        str(_POC_DIR / "data" / "NH3_terminal_His_proton_updated.pdb"),
    )
    ASSEMBLED_PDB = Path(
        environ.get(
            "CONJUGATION_ASSEMBLED_PDB", str(_POC_DIR / "output" / "conjugate_assembled.pdb")
        )
    )
    MINIMIZED_PDB = Path(
        environ.get(
            "CONJUGATION_MINIMIZED_PDB", str(_POC_DIR / "output" / "conjugate_minimized.pdb")
        )
    )
    EQUILIBRATED_PDB = Path(
        environ.get(
            "CONJUGATION_EQUILIBRATED_PDB",
            str(_POC_DIR / "output" / "conjugate_equilibrated.pdb"),
        )
    )
    POLYMER_CACHE = Path(environ.get("CONJUGATION_POLYMER_CACHE", str(_POC_DIR / ".polymer_cache")))
    MONOMER_GROUP_JSON = POLYMER_CACHE / "SBMA-EGPMA-NHS_monomer_group.json"

    _DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data"
    ATRP_INITIATION_RXN = _DATA_DIR / "reactions" / "atrp_initiation.rxn"
    ATRP_POLYMERIZATION_RXN = _DATA_DIR / "reactions" / "atrp_polymerization.rxn"
    ATRP_TERMINATION_RXN = _DATA_DIR / "reactions" / "atrp_termination.rxn"

    # ── Monomer definitions ────────────────────────────────────────────────
    MONOMER_SMILES = {
        "SBMA": "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]",
        "EGPMA": "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]",
        "NHS": "CC(=C)C(=O)ON1C(=O)CCC1=O",
    }
    MONOMER_NAMES = {"A": "SBMA", "B": "EGPMA", "C": "NHS"}
    MONOMER_PROBABILITIES = {"A": 0.945, "B": 0.045, "C": 0.01}
    RESIDUE_NAMES = {"SBMA": "SB", "EGPMA": "EG", "NHS": "NH"}
    NHS_LABEL = "C"

    # ── Conjugation parameters ─────────────────────────────────────────────
    CONJUGATION_SITE_RESIDS = (23, 44)
    N_CONJUGATED = 2
    CONJUGATED_LENGTH = 10
    CENTER_INDEX = CONJUGATED_LENGTH // 2

    # ── Free polymer parameters ────────────────────────────────────────────
    N_FREE_POLYMERS = 5
    FREE_POLYMER_LENGTH = 5

    # ── Placement parameters ───────────────────────────────────────────────
    PACKMOL_TOLERANCE = 2.0
    PACKMOL_REACTIVE_SPHERE_RADIUS = 5.0
    PACKMOL_NLOOP = 500
    PACKMOL_MOVEBADRANDOM = True
    COLLISION_DISTANCE = 1.5
    TARGET_BOND_LENGTH = 1.33
    TARGET_BOND_ANGLE = 120.0
    RNG_SEED = 42

    # ── Free polymer placement ─────────────────────────────────────────────
    FREE_SHELL_MIN_A = 15.0
    FREE_SHELL_MAX_A = 35.0
    N_FREE_PLACEMENT_TRIES = 200
    FREE_COLLISION_DISTANCE = 2.0

    # ── Force field / minimization ─────────────────────────────────────────
    PROTEIN_FF = "ff14sb_off_impropers_0.0.4.offxml"
    SMALL_MOL_FF = "openff-2.0.0.offxml"
    MINIMIZATION_MAX_ITER = 500
    LINKAGE_NEIGHBORHOOD_BONDS = 3

    # ── Vacuum equilibration parameters ────────────────────────────────────
    PROTEIN_HEAVY_RESTRAINT_K = 5000.0  # kJ/mol/nm^2
    VACUUM_EQ_TEMP_K = 310.0
    VACUUM_EQ_TIMESTEP_FS = 1.0
    VACUUM_EQ_STEPS = 10000  # 10 ps
    VACUUM_EQ_FRICTION_PER_PS = 1.0
    PROTEIN_RMSD_THRESHOLD_A = 1.0

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

    @dataclass(frozen=True)
    class FreePlacementResult:
        """Record of a successfully placed free polymer.

        Attributes
        ----------
        polymer_id : str
            Identifier of the free polymer
        sequence : str
            Monomer sequence used for this polymer
        off_molecule : Molecule
            OpenFF molecule carrying the placed conformer
        placed_coords : np.ndarray
            Final placed Cartesian coordinates in angstrom
        """

        polymer_id: str
        sequence: str
        off_molecule: Molecule
        placed_coords: np.ndarray

    def generate_weighted_sequence(length: int, rng: np.random.Generator) -> str:
        """Generate a random polymer sequence from monomer probabilities.

        Parameters
        ----------
        length : int
            Number of residues to sample
        rng : np.random.Generator
            Random generator used for weighted sampling

        Returns
        -------
        str
            Sequence string in monomer label alphabet
        """

        labels = list(MONOMER_PROBABILITIES.keys())
        probs = [MONOMER_PROBABILITIES[x] for x in labels]
        return "".join(rng.choice(labels, size=length, p=probs))

    def generate_conjugated_sequence(length: int, rng: np.random.Generator) -> str:
        """Generate a conjugated polymer sequence with forced center NHS.

        Parameters
        ----------
        length : int
            Number of residues to sample
        rng : np.random.Generator
            Random generator used for weighted sampling

        Returns
        -------
        str
            Sequence string where the center residue is NHS
        """

        seq = list(generate_weighted_sequence(length, rng))
        seq[CENTER_INDEX] = NHS_LABEL
        return "".join(seq)

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
        ATRP_INITIATION_RXN,
        ATRP_POLYMERIZATION_RXN,
        ATRP_TERMINATION_RXN,
        CENTER_INDEX,
        Chem,
        COLLISION_DISTANCE,
        CONJUGATED_LENGTH,
        CONJUGATION_SITE_RESIDS,
        EQUILIBRATED_PDB,
        ForceField,
        FREE_COLLISION_DISTANCE,
        FREE_POLYMER_LENGTH,
        FREE_SHELL_MAX_A,
        FREE_SHELL_MIN_A,
        FragmentGenerator,
        FreePlacementResult,
        KDTree,
        LINKAGE_NEIGHBORHOOD_BONDS,
        MINIMIZATION_MAX_ITER,
        MINIMIZED_PDB,
        MONOMER_GROUP_JSON,
        MONOMER_NAMES,
        MONOMER_PROBABILITIES,
        MONOMER_SMILES,
        MonomerGroup,
        Molecule,
        N_CONJUGATED,
        N_FREE_PLACEMENT_TRIES,
        N_FREE_POLYMERS,
        NHS_LABEL,
        PACKMOL_MOVEBADRANDOM,
        PACKMOL_NLOOP,
        PACKMOL_REACTIVE_SPHERE_RADIUS,
        PACKMOL_TOLERANCE,
        POLYMER_CACHE,
        PROTEIN_FF,
        PROTEIN_HEAVY_RESTRAINT_K,
        PROTEIN_PDB,
        PROTEIN_RMSD_THRESHOLD_A,
        Path,
        PlacementResult,
        PolymerGenerator,
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
        VACUUM_EQ_FRICTION_PER_PS,
        VACUUM_EQ_STEPS,
        VACUUM_EQ_TEMP_K,
        VACUUM_EQ_TIMESTEP_FS,
        generate_conjugated_sequence,
        generate_weighted_sequence,
        logger,
        np,
        unit,
        write_conect_line,
        write_pdb_line,
    )


@app.cell
def _(
    CONJUGATION_SITE_RESIDS,
    ForceField,
    PROTEIN_FF,
    PROTEIN_PDB,
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
    for _resid in CONJUGATION_SITE_RESIDS:
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

    # Suppress per-atom INFO logging from OpenFF's LibraryCharges handler.
    # It emits one line per atom during charge assignment, flooding marimo output.
    import logging as _logging

    _nonbonded_logger = _logging.getLogger("openff.interchange.smirnoff._nonbonded")
    _prev_level = _nonbonded_logger.level
    _nonbonded_logger.setLevel(_logging.WARNING)
    try:
        protein_interchange = protein_ff.create_interchange(protein_top)
    finally:
        _nonbonded_logger.setLevel(_prev_level)
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
        reactive_sites,
    )


@app.cell
def _(
    ATRP_INITIATION_RXN,
    ATRP_POLYMERIZATION_RXN,
    ATRP_TERMINATION_RXN,
    CENTER_INDEX,
    Chem,
    CONJUGATED_LENGTH,
    CONJUGATION_SITE_RESIDS,
    FREE_POLYMER_LENGTH,
    FragmentGenerator,
    MONOMER_GROUP_JSON,
    MONOMER_NAMES,
    MONOMER_SMILES,
    MonomerGroup,
    Molecule,
    N_FREE_POLYMERS,
    NHS_LABEL,
    POLYMER_CACHE,
    PolymerGenerator,
    PolymerPrototype,
    RNG_SEED,
    RESIDUE_NAMES,
    generate_conjugated_sequence,
    generate_weighted_sequence,
    logger,
    np,
    unit,
):
    """Cell 3: Build conjugated/free polymers and charge conjugated prototypes."""

    if not MONOMER_GROUP_JSON.exists():
        logger.info("Generating 3-monomer MonomerGroup (SBMA + EGPMA + NHS)...")
        fragment_gen = FragmentGenerator(
            initiation_rxn_path=ATRP_INITIATION_RXN,
            polymerization_rxn_path=ATRP_POLYMERIZATION_RXN,
            termination_rxn_path=ATRP_TERMINATION_RXN,
            cache_directory=POLYMER_CACHE,
        )
        _monomer_group = fragment_gen.load_or_generate(
            monomer_smiles=MONOMER_SMILES,
            type_prefix="SBMA-EGPMA-NHS",
        )
    else:
        _monomer_group = MonomerGroup.from_file(MONOMER_GROUP_JSON)

    _polymer_gen = PolymerGenerator(
        monomer_group=_monomer_group,
        cache_directory=POLYMER_CACHE,
    )

    rng = np.random.default_rng(RNG_SEED)

    conjugated_sequences = []
    for resid in CONJUGATION_SITE_RESIDS:
        seq = generate_conjugated_sequence(CONJUGATED_LENGTH, rng)
        conjugated_sequences.append(seq)
        logger.info("Conjugated polymer for LYS %d: sequence=%s", resid, seq)

    free_sequences = []
    for i in range(N_FREE_POLYMERS):
        seq = generate_weighted_sequence(FREE_POLYMER_LENGTH, rng)
        free_sequences.append(seq)
        logger.info("Free polymer %d: sequence=%s", i, seq)

    cap_charge_map = {}
    conjugated_prototypes = []

    for _ci, (seq, _resid) in enumerate(
        zip(conjugated_sequences, CONJUGATION_SITE_RESIDS, strict=False)
    ):
        _ = _polymer_gen.generate_polymer(
            sequence=seq,
            monomer_names=MONOMER_NAMES,
            residue_names=RESIDUE_NAMES,
        )

        _polymer_filename = _polymer_gen._make_polymer_filename(seq, MONOMER_NAMES, charged=False)
        polymer_pdb_path = POLYMER_CACHE / f"{_polymer_filename}.pdb"
        polymer_sdf_path = polymer_pdb_path.with_suffix(".sdf")
        if not polymer_pdb_path.exists() or not polymer_sdf_path.exists():
            raise FileNotFoundError(
                f"Missing cached polymer files for sequence {seq}: "
                f"{polymer_pdb_path} / {polymer_sdf_path}"
            )

        logger.info("Conjugated polymer %d (LYS %d): loading %s", _ci, _resid, polymer_pdb_path)

        polymer_rdmol = Chem.MolFromPDBFile(str(polymer_pdb_path), removeHs=False, sanitize=True)
        if polymer_rdmol is None:
            raise ValueError(f"RDKit failed to load {polymer_pdb_path}")

        _suppl = Chem.SDMolSupplier(str(polymer_sdf_path), removeHs=False, sanitize=False)
        _polymer_bo = max((m for m in _suppl if m is not None), key=lambda m: m.GetNumAtoms())
        Chem.SanitizeMol(_polymer_bo)

        polymer_off = Molecule.from_rdkit(
            polymer_rdmol,
            allow_undefined_stereo=True,
            hydrogens_are_explicit=True,
        )

        # Middle residues start at PDB residue 1 because polymerist appends
        # terminal 1-site fragments at the end, so sequence index 5 maps to residue 5
        reactive_residue_number = CENTER_INDEX
        _nh2_indices = set()
        for _atom in polymer_rdmol.GetAtoms():
            _info = _atom.GetPDBResidueInfo()
            if (
                _info
                and _info.GetResidueName().strip().startswith(RESIDUE_NAMES["NHS"])
                and _info.GetResidueNumber() == reactive_residue_number
            ):
                _nh2_indices.add(_atom.GetIdx())

        if not _nh2_indices:
            raise ValueError(
                f"No NHS residue at center (residue {reactive_residue_number}) for sequence {seq}"
            )

        _reactive_c = None
        _ester_o = None
        for _idx in _nh2_indices:
            _atom = polymer_rdmol.GetAtomWithIdx(_idx)
            if _atom.GetSymbol() != "C":
                continue
            _o_nbrs = [n.GetIdx() for n in _atom.GetNeighbors() if n.GetSymbol() == "O"]
            if len(_o_nbrs) == 2:
                _reactive_c = _idx
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
        _retained_sorted = tuple(i for i in range(polymer_rdmol.GetNumAtoms()) if i not in _leaving)

        _core_rw = Chem.RWMol(_polymer_bo)
        for _idx in sorted(_leaving, reverse=True):
            _core_rw.RemoveAtom(_idx)
        _core = _core_rw.GetMol()
        Chem.SanitizeMol(_core)

        _old_to_core = {}
        _cursor = 0
        for _old_idx in range(polymer_rdmol.GetNumAtoms()):
            if _old_idx not in _leaving:
                _old_to_core[_old_idx] = _cursor
                _cursor += 1
        _rc_core = _old_to_core[_reactive_c]

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

        _merged_for_charge = Chem.CombineMols(_core, _cap_rd)
        _merged_rw = Chem.RWMol(_merged_for_charge)
        _cap_offset = _core.GetNumAtoms()
        _merged_rw.AddBond(_rc_core, _cap_offset, Chem.BondType.SINGLE)
        _merged_rw.GetAtomWithIdx(_cap_offset).SetNumRadicalElectrons(0)
        _product_rd = _merged_rw.GetMol()
        Chem.SanitizeMol(_product_rd)

        _product_off = Molecule.from_rdkit(
            _product_rd,
            allow_undefined_stereo=True,
            hydrogens_are_explicit=True,
        )

        try:
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

            _product_off.assign_partial_charges(
                partial_charge_method="openff-gnn-am1bcc-0.1.0-rc.3.pt",
                toolkit_registry=NAGLToolkitWrapper(),
            )
        except (ImportError, ValueError, RuntimeError) as _exc:
            logger.warning("NAGL failed (%s), falling back to am1bcc", _exc)
            _product_off.assign_partial_charges(partial_charge_method="am1bcc")

        _product_charges_all = np.array(
            [c.m_as(unit.elementary_charge) for c in _product_off.partial_charges]
        )
        _product_charges_core = _product_charges_all[: _core.GetNumAtoms()]
        _retained_charges = np.array(
            [_product_charges_core[_old_to_core[i]] for i in _retained_sorted]
        )

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

        _cap_n_atom = _product_off.atom(_cap_offset)
        for _bond in _cap_n_atom.bonds:
            _nbr_idx = _bond.atom1_index if _bond.atom2_index == _cap_offset else _bond.atom2_index
            if _product_off.atom(_nbr_idx).atomic_number == 1:
                _cap_charge_map["HZ1"] = float(_product_charges_all[_nbr_idx])
                break

        if not cap_charge_map:
            cap_charge_map = _cap_charge_map

        conjugated_prototypes.append(
            PolymerPrototype(
                off_molecule=polymer_off,
                rdmol=polymer_rdmol,
                reactive_carbon_idx=_reactive_c,
                leaving_group_indices=_leaving_sorted,
                retained_atom_indices=_retained_sorted,
                product_charges=_retained_charges,
            )
        )

    free_polymer_offs = []
    for _fi, seq in enumerate(free_sequences):
        charged_off = _polymer_gen.generate_polymer(
            sequence=seq,
            monomer_names=MONOMER_NAMES,
            residue_names=RESIDUE_NAMES,
        )
        _free_off = Molecule(charged_off)
        free_polymer_offs.append(_free_off)
        logger.info(
            "Free polymer %d: %d atoms, charge=%.4f",
            _fi,
            _free_off.n_atoms,
            sum(c.m_as(unit.elementary_charge) for c in _free_off.partial_charges),
        )

    return (
        cap_charge_map,
        conjugated_prototypes,
        conjugated_sequences,
        free_polymer_offs,
        free_sequences,
    )


@app.cell
def _(
    CONJUGATION_SITE_RESIDS,
    conjugated_prototypes,
    logger,
    reactive_sites: "dict[tuple[str, int], ReactiveSite]",
):
    """Cell 4: Assign conjugated polymer prototypes to fixed lysine sites."""

    assignments = []
    for _ci, _resid in enumerate(CONJUGATION_SITE_RESIDS):
        site = reactive_sites[("A", _resid)]
        prototype = conjugated_prototypes[_ci]
        assignments.append((site, prototype))
        logger.info(
            "  LYS %d <- conjugated polymer %d (%d atoms)",
            _resid,
            _ci,
            prototype.off_molecule.n_atoms,
        )
    return (assignments,)


@app.cell
def _(
    FREE_COLLISION_DISTANCE,
    FREE_SHELL_MAX_A,
    FREE_SHELL_MIN_A,
    FreePlacementResult,
    KDTree,
    Molecule,
    N_FREE_PLACEMENT_TRIES,
    PACKMOL_MOVEBADRANDOM,
    PACKMOL_NLOOP,
    PACKMOL_REACTIVE_SPHERE_RADIUS,
    PACKMOL_TOLERANCE,
    Path,
    PlacementResult,
    RNG_SEED,
    Rotation,
    TARGET_BOND_LENGTH,
    assignments: list[tuple],
    centroid_vec,
    free_polymer_offs,
    free_sequences,
    logger,
    np,
    protein_coords,
    reactive_sites: "dict[tuple[str, int], ReactiveSite]",
):
    """Cell 5: Packmol-based polymer placement at lysine conjugation sites.

    Uses Packmol to place each conjugated polymer near its target lysine NZ,
    with per-atom sphere constraints on the reactive carbon. After Packmol
    packing, a post-placement translation snaps the reactive carbon to
    TARGET_BOND_LENGTH from NZ.

    Free (non-covalent) polymers use shell-based random placement (unchanged).
    """
    import tempfile

    from polyzymd.utils.packmol import build_packmol_input, run_packmol

    # ── Helper: write a PDB from coordinates + elements ────────────────────
    def _write_simple_pdb(path: Path, coords: np.ndarray, elements: list[str]) -> None:
        """Write a minimal PDB file from coordinates and element symbols."""
        lines = []
        for i, (xyz, elem) in enumerate(zip(coords, elements, strict=False)):
            name = f"{elem}{i + 1:>3d}"[:4]
            lines.append(
                f"HETATM{i + 1:>5d} {name:<4s} UNK A   1    "
                f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}"
                f"  1.00  0.00          {elem:>2s}\n"
            )
        lines.append("END\n")
        path.write_text("".join(lines))

    # ── Helper: Kabsch alignment ───────────────────────────────────────────
    def _kabsch_align(P: np.ndarray, Q: np.ndarray):
        """Compute rotation R and translation t that minimizes ||R@P + t - Q||.

        Parameters
        ----------
        P : (N, 3) source points (template)
        Q : (N, 3) target points (packed)

        Returns
        -------
        R : (3, 3) rotation matrix
        t : (3,) translation vector
        """
        centroid_p = P.mean(axis=0)
        centroid_q = Q.mean(axis=0)
        P_c = P - centroid_p
        Q_c = Q - centroid_q
        H = P_c.T @ Q_c
        U, _, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        S = np.diag([1.0, 1.0, d])  # correct for reflection
        R = Vt.T @ S @ U.T
        t = centroid_q - R @ centroid_p
        return R, t

    # ── Prepare protein PDB for Packmol (exclude NZ + removable HZ) ───────
    # Exclude NZ and the 2 HZ atoms being removed at each conjugation site
    # This prevents Packmol from over-blocking the bond-forming region
    exclude_indices: set[int] = set()
    for _site_key in reactive_sites:
        _site_obj = reactive_sites[_site_key]
        exclude_indices.add(_site_obj.nz_index)
        for _hz in _site_obj.hz_indices[1:]:
            exclude_indices.add(_hz)

    keep_mask = np.ones(len(protein_coords), dtype=bool)
    keep_mask[list(exclude_indices)] = False
    trimmed_protein_coords = protein_coords[keep_mask]

    # Get element symbols for trimmed protein (needed for PDB writing)
    # We need to recover the element for each kept atom. Since protein_coords
    # come from protein_mol (Cell 2), we can use the same OpenFF molecule
    # But protein_mol isn't passed to Cell 5. Instead, we use a simple
    # heuristic: write all as carbon (Packmol only needs coords for steric
    # exclusion; element doesn't affect packing)
    trimmed_elements = ["C"] * len(trimmed_protein_coords)

    # ── Prepare polymer PDBs (retained atoms only, centered on reactive C) ─
    polymer_templates = []  # list of dicts with template info
    for site, polymer in assignments:
        retained_idx = np.asarray(polymer.retained_atom_indices, dtype=int)
        # Use the cached SDF conformer (from polymer generation)
        full_coords = polymer.off_molecule.conformers[0].m_as("angstrom")
        retained_coords = full_coords[retained_idx]
        reactive_c_local = np.where(retained_idx == polymer.reactive_carbon_idx)[0][0]

        # Elements for retained atoms
        elements = [polymer.off_molecule.atoms[i].symbol for i in retained_idx]

        polymer_templates.append(
            {
                "site": site,
                "polymer": polymer,
                "full_coords": full_coords,
                "retained_idx": retained_idx,
                "retained_coords": retained_coords,
                "reactive_c_local": reactive_c_local,
                "elements": elements,
            }
        )

    # ── Compute bounding box + coordinate shift ────────────────────────────
    # Packmol requires positive coordinates. Compute bounding box of all
    # components and shift everything to the positive octant
    all_coords = [trimmed_protein_coords]
    for tmpl in polymer_templates:
        all_coords.append(tmpl["retained_coords"])
    combined = np.vstack(all_coords)
    global_min = combined.min(axis=0)
    # Add padding so nothing sits at exactly 0
    SHIFT_PADDING = 10.0
    coord_shift = -global_min + SHIFT_PADDING

    # Shifted protein coords
    shifted_protein = trimmed_protein_coords + coord_shift

    # Shifted NZ positions (for sphere constraints)
    shifted_nz = {}
    for site, _polymer in assignments:
        shifted_nz[site.resid] = site.nz_position + coord_shift

    # Box size: bounding box of shifted coords + generous padding
    BOX_PADDING = 30.0
    shifted_all = np.vstack(
        [shifted_protein] + [tmpl["retained_coords"] + coord_shift for tmpl in polymer_templates]
    )
    box_size = shifted_all.max(axis=0) + BOX_PADDING

    # ── Write PDB files to temp directory ──────────────────────────────────
    with tempfile.TemporaryDirectory(prefix="packmol_conjugation_") as _tmp_dir:
        work_dir = Path(_tmp_dir)
        logger.info("Packmol working directory: %s", work_dir)

        # Write trimmed protein PDB
        protein_pdb = work_dir / "protein_trimmed.pdb"
        _write_simple_pdb(protein_pdb, shifted_protein, trimmed_elements)

        # Write polymer PDBs (shifted, NOT centered on reactive C — Packmol
        # handles placement. We keep the original relative geometry)
        polymer_pdbs = []
        structure_extra_lines = []
        for i, tmpl in enumerate(polymer_templates):
            shifted_retained = tmpl["retained_coords"] + coord_shift
            pdb_path = work_dir / f"polymer_{i}.pdb"
            _write_simple_pdb(pdb_path, shifted_retained, tmpl["elements"])
            polymer_pdbs.append(str(pdb_path))

            # Per-atom constraint: reactive carbon inside sphere near NZ
            nz_shifted = shifted_nz[tmpl["site"].resid]
            # Packmol uses 1-based indexing
            rc_1based = tmpl["reactive_c_local"] + 1
            extra = [
                f"atoms {rc_1based}",
                (
                    "inside sphere "
                    f"{nz_shifted[0]:.6f} {nz_shifted[1]:.6f} {nz_shifted[2]:.6f} "
                    f"{PACKMOL_REACTIVE_SPHERE_RADIUS:.1f}"
                ),
                "end atoms",
            ]
            structure_extra_lines.append(extra)

        # ── Build and run Packmol ──────────────────────────────────────────
        input_text = build_packmol_input(
            molecule_pdb_paths=polymer_pdbs,
            molecule_counts=[1] * len(polymer_pdbs),
            box_size_angstrom=box_size,
            tolerance_angstrom=PACKMOL_TOLERANCE,
            solute_pdb_path=str(protein_pdb),
            use_pbc=False,
            movebadrandom=PACKMOL_MOVEBADRANDOM,
            nloop=PACKMOL_NLOOP,
            structure_extra_lines=structure_extra_lines,
        )

        logger.info("Packmol input:\n%s", input_text)
        output_pdb = run_packmol(input_text, work_dir)
        logger.info("Packmol output: %s", output_pdb)

        # ── Parse Packmol output ───────────────────────────────────────────
        # Read all ATOM/HETATM lines from the output PDB
        packed_coords_all = []
        with open(output_pdb, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    packed_coords_all.append([x, y, z])
        packed_coords_all = np.array(packed_coords_all)

        n_protein = len(shifted_protein)
        expected_packed_atoms = n_protein + sum(
            len(tmpl["retained_coords"]) for tmpl in polymer_templates
        )
        if len(packed_coords_all) != expected_packed_atoms:
            raise RuntimeError(
                f"Packmol output has {len(packed_coords_all)} atoms, "
                f"expected {expected_packed_atoms} (protein={n_protein} + "
                f"polymers={expected_packed_atoms - n_protein})"
            )

        # Slice: first N_protein atoms are the fixed protein, then each polymer
        polymer_packed_coords = []
        offset = n_protein
        for tmpl in polymer_templates:
            n_atoms = len(tmpl["retained_coords"])
            polymer_packed_coords.append(packed_coords_all[offset : offset + n_atoms])
            offset += n_atoms

    # ── Kabsch alignment + post-placement refinement ───────────────────────
    # For each polymer:
    # 1. Kabsch-align the original retained template → packed retained coords
    # 2. Apply the same transform to the FULL original coords
    # 3. Translate so reactive C is at TARGET_BOND_LENGTH from NZ
    # 4. Undo the coordinate shift

    occupied = trimmed_protein_coords.copy()  # back in original frame
    placement_results: list[PlacementResult] = []

    for i, tmpl in enumerate(polymer_templates):
        # Source: original retained coords (before shift)
        P = tmpl["retained_coords"]
        # Target: packed coords (undo shift)
        Q = polymer_packed_coords[i] - coord_shift

        R, t = _kabsch_align(P, Q)

        # Apply transform to FULL original coords
        full_transformed = (R @ tmpl["full_coords"].T).T + t

        # Post-placement: translate so reactive C is exactly at
        # TARGET_BOND_LENGTH from NZ
        nz_pos = tmpl["site"].nz_position
        rc_idx = tmpl["polymer"].reactive_carbon_idx
        current_rc = full_transformed[rc_idx]
        nz_to_rc = current_rc - nz_pos
        nz_to_rc_dist = np.linalg.norm(nz_to_rc)
        if nz_to_rc_dist < 1e-8:
            raise RuntimeError(f"Reactive carbon coincides with NZ for LYS {tmpl['site'].resid}")
        nz_to_rc_hat = nz_to_rc / nz_to_rc_dist
        target_rc = nz_pos + TARGET_BOND_LENGTH * nz_to_rc_hat
        translation = target_rc - current_rc
        full_transformed += translation

        # Informational sanity check after snapping to bond length
        protein_tree = KDTree(protein_coords)
        snap_dists, _ = protein_tree.query(full_transformed)
        min_post_snap_protein_dist = float(np.min(snap_dists))
        if min_post_snap_protein_dist < 1.0:
            logger.warning(
                "LYS %d: post-snap minimum polymer-protein distance %.3f A (< 1.0 A)",
                tmpl["site"].resid,
                min_post_snap_protein_dist,
            )

        # Check clearance against occupied atoms
        retained_idx = tmpl["retained_idx"]
        collision_mask = retained_idx != tmpl["polymer"].reactive_carbon_idx
        check_pts = full_transformed[retained_idx[collision_mask]]
        if len(occupied) > 0:
            tree = KDTree(occupied)
            dists, _ = tree.query(check_pts)
            min_dist = float(np.min(dists))
        else:
            min_dist = float("inf")

        logger.info(
            "LYS %d: Packmol placed, post-snap distance to NZ=%.3f A, min clearance=%.2f A",
            tmpl["site"].resid,
            np.linalg.norm(full_transformed[rc_idx] - nz_pos),
            min_dist,
        )

        placement_results.append(
            PlacementResult(
                site=tmpl["site"],
                polymer=tmpl["polymer"],
                placed_coords=full_transformed,
                min_protein_distance=min_dist,
            )
        )

        # Add retained coords to occupied set for next polymer
        occupied = np.vstack((occupied, full_transformed[retained_idx]))

    # ── Free polymer placement (unchanged from original) ───────────────────
    free_placement_results: list[FreePlacementResult] = []
    for _fi, _free_off in enumerate(free_polymer_offs):
        rng_free = np.random.default_rng(RNG_SEED + 1000 + _fi)

        if not _free_off.conformers:
            raise RuntimeError(f"No conformers for free polymer {_fi}")

        placed = False
        for _ in range(N_FREE_PLACEMENT_TRIES):
            conf_coords = _free_off.conformers[0].m_as("angstrom")

            centered = conf_coords - conf_coords.mean(axis=0)
            rot = Rotation.random(random_state=rng_free.integers(0, 2**31))
            rotated = rot.apply(centered)

            direction = rng_free.standard_normal(3)
            direction /= np.linalg.norm(direction)
            radius = rng_free.uniform(FREE_SHELL_MIN_A, FREE_SHELL_MAX_A)
            _offset = centroid_vec + radius * direction
            trial = rotated + _offset

            occupied_tree = KDTree(occupied)
            dists, _ = occupied_tree.query(trial)
            min_dist = float(np.min(dists))
            if min_dist >= FREE_COLLISION_DISTANCE:
                from openff.units import Quantity as Q
                from openff.units import unit as u

                _free_off._conformers = [Q(trial, u.angstrom)]

                occupied = np.vstack((occupied, trial))
                free_placement_results.append(
                    FreePlacementResult(
                        polymer_id=f"free_{_fi}",
                        sequence=free_sequences[_fi],
                        off_molecule=_free_off,
                        placed_coords=trial,
                    )
                )
                placed = True
                logger.info("Placed free polymer %d: clearance %.2f A", _fi, min_dist)
                break

        if not placed:
            raise RuntimeError(
                f"Failed to place free polymer {_fi} after {N_FREE_PLACEMENT_TRIES} tries"
            )

    placed_free_polymer_offs = [res.off_molecule for res in free_placement_results]

    return free_placement_results, placed_free_polymer_offs, placement_results


@app.cell
def _(
    Chem,
    Molecule,
    PROTEIN_PDB,
    Quantity,
    cap_charge_map,
    logger,
    np,
    placed_free_polymer_offs,
    placement_results: "list[PlacementResult]",
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
        PROTEIN_PDB,
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
    _n_prot_orig = _prot_rd.GetNumAtoms()

    _prot_removed: set[int] = set()  # original protein indices removed so far

    for _pi, _placement in enumerate(placement_results):
        _site = _placement.site
        _polymer = _placement.polymer

        # ── Identify which 2 HZ to remove from this lysine ────────────────
        # Keep the first HZ, remove the other two.
        # hz_indices are in the ORIGINAL protein_mol ordering.
        _hz_to_remove = list(_site.hz_indices[1:])  # remove last 2, keep first
        _prot_removed.update(_hz_to_remove)

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
    # Build protein charge array (remove HZ entries, redistribute charge)
    _prot_charges_mod: list[float] = []
    for _old in range(_n_prot_orig):
        if _old not in _prot_removed:
            _prot_charges_mod.append(float(protein_charge_arr[_old]))

    # Build protein coordinate array (remove H entries)
    _prot_coords_mod = protein_coords[
        np.array([_old not in _prot_removed for _old in range(_n_prot_orig)], dtype=bool)
    ]

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
            _pi,
            _site.resid,
            _merged.GetNumAtoms(),
        )

    # ── Patch lysine junction charges ──────────────────────────────────────
    # For each conjugated lysine, replace the ff14SB charges on
    # NZ, CE, CD, CG, CB with the cap-derived charges.  Also patch the
    # retained NZ hydrogen (HZ1 equivalent) — identified by index rather
    # than name, since OpenFF/CCD may use non-standard H names (e.g. H10).
    for _placement in placement_results:
        _site = _placement.site
        for _atom_name, _cap_q in cap_charge_map.items():
            if _atom_name == "HZ1":
                # The retained HZ is the first in hz_indices (others removed).
                _hz_orig_idx = _site.hz_indices[0]
                if _hz_orig_idx in _prot_old_to_new:
                    _new_idx = _prot_old_to_new[_hz_orig_idx]
                    _old_q = _merged_charges[_new_idx]
                    _merged_charges[_new_idx] = _cap_q
                    logger.info(
                        "  Patched LYS %d HZ1 (atom %d): %.4f → %.4f",
                        _site.resid,
                        _hz_orig_idx,
                        _old_q,
                        _cap_q,
                    )
                continue
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
                            _site.resid,
                            _atom_name,
                            _old_q,
                            _cap_q,
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
        _merged_noconf,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )
    conjugate_off.add_conformer(Quantity(_merged_coords, unit.angstrom))

    # ── Fix formal charges on protein atoms ────────────────────────────────
    # RDKit assigns different formal charges than OpenFF/CCD for some residues
    # (e.g. ARG NH1 gets FC=+1 in RDKit but FC=0 in OpenFF, and HIS ND1 gets
    # FC=+1 for aromatic perception while CCD uses FC=0).
    #
    # We correct most of them, but SKIP aromatic atoms (e.g. HIS ND1) because
    # create_interchange() roundtrips through RDKit, which would reject the
    # change as a valence violation.  The construction-accounting charge target
    # below accounts for these uncorrectable FC mismatches.
    #
    # Conjugated lysine NZ changes from +1 (protonated NH3+) to
    # 0 (amide N) after removing 2 HZ and forming the crosslink bond.
    _conjugated_nz_indices = {_prot_old_to_new[p.site.nz_index] for p in placement_results}
    _n_fc_fixed = 0
    _fc_skipped: list[tuple[int, int, int]] = []  # (new_idx, current_fc, target_fc)
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
        # Safety: don't change FC on aromatic atoms — create_interchange()
        # roundtrips through RDKit which enforces kekulization/valence rules.
        _rd_atom = _merged.GetAtomWithIdx(_new_idx)
        if _rd_atom.GetIsAromatic():
            _fc_skipped.append((_new_idx, _current_fc, _target_fc))
            logger.info(
                "  Skipping FC fix for aromatic atom %d (%s %s %s): "
                "keeping FC=%d (target=%d) — will absorb in partial charges",
                _orig_idx,
                protein_mol.atom(_orig_idx).metadata.get("residue_name", "?"),
                protein_mol.atom(_orig_idx).metadata.get("residue_number", "?"),
                protein_mol.atom(_orig_idx).name,
                _current_fc,
                _target_fc,
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
    logger.info("Fixed %d formal charges (%d skipped aromatic)", _n_fc_fixed, len(_fc_skipped))

    # ── Compute target net charge from construction accounting ────────────
    # We derive the intended system charge from first principles rather than
    # trusting sum(formal_charge) on the merged molecule, as a cross-check.
    _protein_formal_charge = sum(a.formal_charge.m for a in protein_mol.atoms)
    _conjugation_fc_delta = -1 * len(placement_results)
    _polymer_fc_total = 0
    for _placement in placement_results:
        _poly_off = _placement.polymer.off_molecule
        _retained = set(_placement.polymer.retained_atom_indices)
        _polymer_fc_total += sum(_poly_off.atom(i).formal_charge.m for i in _retained)
    _target_charge = _protein_formal_charge + _conjugation_fc_delta + _polymer_fc_total

    # The molecule's formal charge (as seen by create_interchange) may differ
    # from _target_charge if aromatic atoms had uncorrectable FC mismatches.
    # create_interchange requires partial charges to sum to the molecule FC,
    # so we must target that value.  But we handle the discrepancy surgically:
    # for each skipped aromatic atom, absorb the FC delta directly into its
    # partial charge (rather than smearing it across all atoms).
    _molecule_fc = sum(a.formal_charge.m for a in conjugate_off.atoms)
    _fc_discrepancy = _molecule_fc - _target_charge

    # ── Charge ledger ──────────────────────────────────────────────────────
    _charge_arr = np.array(_merged_charges)
    _charge_sum = _charge_arr.sum()
    logger.info("── Charge ledger ──")
    logger.info("  Protein formal charge (from CCD):  %+.0f", _protein_formal_charge)
    logger.info(
        "  Conjugation FC delta (%d sites):   %+.0f",
        len(placement_results),
        _conjugation_fc_delta,
    )
    logger.info(
        "  Polymer formal charge (%d chains):  %+.0f",
        len(placement_results),
        _polymer_fc_total,
    )
    logger.info("  Target net charge (physical):       %+.0f", _target_charge)
    logger.info("  Molecule FC (after corrections):    %+.0f", _molecule_fc)
    if _fc_discrepancy != 0:
        logger.info(
            "  FC discrepancy (uncorrectable):     %+.0f (from %d skipped aromatic atoms)",
            _fc_discrepancy,
            len(_fc_skipped),
        )
    logger.info("  Partial charge sum (pre-fix):       %+.4f", _charge_sum)

    # Step 1: For each skipped aromatic atom, adjust its partial charge by
    # (current_fc − target_fc) so the ff14SB charge reflects the correct
    # physical state despite the incorrect formal charge.
    for _new_idx, _current_fc, _target_fc_val in _fc_skipped:
        _delta = _current_fc - _target_fc_val
        _old_q = _charge_arr[_new_idx]
        _charge_arr[_new_idx] += _delta
        logger.info(
            "  Absorbed FC artifact on atom %d: q %.4f → %.4f (delta %+d)",
            _new_idx,
            _old_q,
            _charge_arr[_new_idx],
            _delta,
        )

    # Step 2: Redistribute the remaining small residual (ff14SB/NAGL junction
    # mismatch) evenly across all atoms to hit the molecule FC exactly.
    _charge_sum_after_absorb = _charge_arr.sum()
    _residual = _molecule_fc - _charge_sum_after_absorb
    logger.info("  After absorbing FC artifacts:       %+.4f", _charge_sum_after_absorb)
    logger.info("  Residual to redistribute:           %+.4f", _residual)

    if abs(_residual) > 1.0:
        raise ValueError(
            "Charge mismatch too large to redistribute safely: "
            f"molecule_fc={_molecule_fc:.4f}, sum={_charge_sum_after_absorb:.4f}, "
            f"residual={_residual:.4f}. Check charge construction."
        )
    if abs(_residual) > 0.01:
        logger.info(
            "  Redistributing %.4f e across %d atoms (%.6f e/atom)",
            _residual,
            len(_charge_arr),
            _residual / len(_charge_arr),
        )
        _charge_arr += _residual / len(_charge_arr)
    else:
        logger.info("  Residual within tolerance (<=0.01 e), no redistribution needed")

    conjugate_off.partial_charges = Quantity(_charge_arr, unit.elementary_charge)

    logger.info(
        "Conjugate molecule: %d atoms, physical charge=%+.0f, "
        "molecule FC=%+.0f, partial charge sum=%.4f",
        conjugate_off.n_atoms,
        _target_charge,
        _molecule_fc,
        _charge_arr.sum(),
    )

    prot_removed_indices = _prot_removed

    protein_atom_indices = sorted(_prot_old_to_new.values())
    protein_heavy_indices = [
        _prot_old_to_new[orig]
        for orig in range(_n_prot_orig)
        if orig not in _prot_removed and protein_mol.atom(orig).atomic_number > 1
    ]

    protein_backbone_heavy_indices = []
    for orig in range(_n_prot_orig):
        if orig in _prot_removed:
            continue
        atom = protein_mol.atom(orig)
        if atom.atomic_number > 1 and atom.name in ("CA", "C", "N", "O"):
            protein_backbone_heavy_indices.append(_prot_old_to_new[orig])

    conjugated_polymer_ranges = []
    _offset = _prot_mod.GetNumAtoms()
    for _pi, _placement in enumerate(placement_results):
        n_retained = len(_placement.polymer.retained_atom_indices)
        conj_indices = list(range(_offset, _offset + n_retained))
        conj_heavy = [i for i in conj_indices if conjugate_off.atom(i).atomic_number > 1]
        conjugated_polymer_ranges.append(
            {
                "polymer_id": f"conj_{_pi}",
                "site_resid": _placement.site.resid,
                "atom_indices": conj_indices,
                "heavy_atom_indices": conj_heavy,
            }
        )
        _offset += n_retained

    all_conjugate_heavy = protein_heavy_indices.copy()
    for cr in conjugated_polymer_ranges:
        all_conjugate_heavy.extend(cr["heavy_atom_indices"])
    all_conjugate_heavy.sort()

    free_polymer_ranges = []
    free_offset = conjugate_off.n_atoms
    for _fi, free_off_mol in enumerate(placed_free_polymer_offs):
        n_free = free_off_mol.n_atoms
        free_indices = list(range(free_offset, free_offset + n_free))
        free_heavy = [
            i for i in free_indices if free_off_mol.atom(i - free_offset).atomic_number > 1
        ]
        free_polymer_ranges.append(
            {
                "polymer_id": f"free_{_fi}",
                "atom_indices": free_indices,
                "heavy_atom_indices": free_heavy,
            }
        )
        free_offset += n_free

    all_free_heavy = []
    for fr in free_polymer_ranges:
        all_free_heavy.extend(fr["heavy_atom_indices"])
    all_free_heavy.sort()

    component_metadata = {
        "n_conjugate_atoms": conjugate_off.n_atoms,
        "n_free_polymer_atoms": sum(m.n_atoms for m in placed_free_polymer_offs),
        "protein": {
            "atom_indices": protein_atom_indices,
            "heavy_atom_indices": protein_heavy_indices,
            "backbone_heavy_atom_indices": protein_backbone_heavy_indices,
        },
        "conjugated_polymers": conjugated_polymer_ranges,
        "free_polymers": free_polymer_ranges,
        "restraint_groups": {
            "protein_heavy": protein_heavy_indices,
            "protein_backbone_heavy": protein_backbone_heavy_indices,
            "conjugate_heavy": all_conjugate_heavy,
            "free_polymer_heavy": all_free_heavy,
        },
    }

    logger.info(
        "Component metadata: %d protein heavy, %d backbone heavy, %d conjugate heavy, "
        "%d free heavy",
        len(protein_heavy_indices),
        len(protein_backbone_heavy_indices),
        len(all_conjugate_heavy),
        len(all_free_heavy),
    )

    return component_metadata, conjugate_off, prot_removed_indices


@app.cell
def _(
    ASSEMBLED_PDB,
    PROTEIN_PDB,
    Path,
    free_placement_results,
    logger,
    np,
    placement_results: "list[PlacementResult]",
    prot_removed_indices,
    write_conect_line,
    write_pdb_line,
):
    """Cell 7: Write assembled PDB following PolyzyMD chain convention."""

    _pdb_path = PROTEIN_PDB
    _atom_records = []
    with open(_pdb_path, "r", encoding="utf-8") as _f:
        for _line in _f:
            if _line.startswith(("ATOM", "HETATM")):
                _atom_records.append(
                    {
                        "line": _line,
                        "serial": int(_line[6:11]),
                        "name": _line[12:16].strip(),
                        "resname": _line[17:20].strip(),
                        "chain": _line[21].strip() or "A",
                        "resid": int(_line[22:26]),
                    }
                )

    _serial_lookup = {}
    _existing_serials = set()
    kept_lines = []
    for _orig_idx, _rec in enumerate(_atom_records):
        if _orig_idx in prot_removed_indices:
            continue
        kept_lines.append(_rec["line"])
        _serial_lookup[_orig_idx] = _rec["serial"]
        _existing_serials.add(_rec["serial"])

    _next_serial = max(_existing_serials) + 1

    polymer_lines = []
    conect_map = {}
    chain_c_resid = 1

    for _pi, _placement in enumerate(placement_results):
        _polymer = _placement.polymer

        local_resids = set()
        for _orig_idx in _polymer.retained_atom_indices:
            _rd_atom = _polymer.rdmol.GetAtomWithIdx(_orig_idx)
            _rd_info = _rd_atom.GetPDBResidueInfo()
            if _rd_info:
                local_resids.add(_rd_info.GetResidueNumber())
        local_resids_sorted = sorted(local_resids)
        local_to_global = {r: chain_c_resid + i for i, r in enumerate(local_resids_sorted)}
        chain_c_resid += len(local_resids_sorted)

        _serial_map = {}
        for _ri, _orig_idx in enumerate(_polymer.retained_atom_indices):
            _atom = _polymer.off_molecule.atoms[_orig_idx]
            _x, _y, _z = _placement.placed_coords[_orig_idx]
            _rd_atom = _polymer.rdmol.GetAtomWithIdx(_orig_idx)
            _rd_info = _rd_atom.GetPDBResidueInfo()
            _atom_name = _rd_info.GetName().strip() if _rd_info else f"{_atom.symbol}{_ri:03d}"
            _res_name = _rd_info.GetResidueName().strip() if _rd_info else "POL"
            _local_resid = _rd_info.GetResidueNumber() if _rd_info else 1

            polymer_lines.append(
                write_pdb_line(
                    serial=_next_serial,
                    atom_name=_atom_name,
                    residue_name=_res_name,
                    chain_id="C",
                    residue_number=local_to_global.get(_local_resid, chain_c_resid - 1),
                    x=float(_x),
                    y=float(_y),
                    z=float(_z),
                    element=_atom.symbol,
                    record="HETATM",
                )
            )
            _serial_map[_orig_idx] = _next_serial
            _next_serial += 1

        _retained_set = set(_polymer.retained_atom_indices)
        for _bond in _polymer.off_molecule.bonds:
            _i = _bond.atom1_index
            _j = _bond.atom2_index
            if _i in _retained_set and _j in _retained_set:
                _si = _serial_map[_i]
                _sj = _serial_map[_j]
                conect_map.setdefault(_si, set()).add(_sj)
                conect_map.setdefault(_sj, set()).add(_si)

        _nz_serial = _serial_lookup[_placement.site.nz_index]
        _rc_serial = _serial_map[_polymer.reactive_carbon_idx]
        conect_map.setdefault(_nz_serial, set()).add(_rc_serial)
        conect_map.setdefault(_rc_serial, set()).add(_nz_serial)

        polymer_lines.append("TER\n")

    for _fp in free_placement_results:
        _off_mol = _fp.off_molecule
        _coords = _fp.placed_coords

        _rd_mol = _off_mol.to_rdkit()
        local_resids = set()
        for _ai in range(_off_mol.n_atoms):
            _rd_atom = _rd_mol.GetAtomWithIdx(_ai)
            _rd_info = _rd_atom.GetPDBResidueInfo()
            if _rd_info:
                local_resids.add(_rd_info.GetResidueNumber())
            else:
                local_resids.add(1)
        local_resids_sorted = sorted(local_resids)
        local_to_global = {r: chain_c_resid + i for i, r in enumerate(local_resids_sorted)}
        chain_c_resid += len(local_resids_sorted)

        _serial_map_free = {}
        for _ai in range(_off_mol.n_atoms):
            _atom = _off_mol.atoms[_ai]
            _x, _y, _z = _coords[_ai]
            _rd_atom = _rd_mol.GetAtomWithIdx(_ai)
            _rd_info = _rd_atom.GetPDBResidueInfo()
            _atom_name = _rd_info.GetName().strip() if _rd_info else f"{_atom.symbol}{_ai:03d}"
            _res_name = _rd_info.GetResidueName().strip() if _rd_info else "POL"
            _local_resid = _rd_info.GetResidueNumber() if _rd_info else 1

            polymer_lines.append(
                write_pdb_line(
                    serial=_next_serial,
                    atom_name=_atom_name,
                    residue_name=_res_name,
                    chain_id="C",
                    residue_number=local_to_global.get(_local_resid, chain_c_resid - 1),
                    x=float(_x),
                    y=float(_y),
                    z=float(_z),
                    element=_atom.symbol,
                    record="HETATM",
                )
            )
            _serial_map_free[_ai] = _next_serial
            _next_serial += 1

        for _bond in _off_mol.bonds:
            _si = _serial_map_free[_bond.atom1_index]
            _sj = _serial_map_free[_bond.atom2_index]
            conect_map.setdefault(_si, set()).add(_sj)
            conect_map.setdefault(_sj, set()).add(_si)

        polymer_lines.append("TER\n")

    assembled_path = Path(ASSEMBLED_PDB)
    with open(assembled_path, "w", encoding="utf-8") as _f:
        for _line in kept_lines:
            _f.write(_line)
        _f.write("TER\n")
        for _line in polymer_lines:
            _f.write(_line)
        for _serial in sorted(conect_map):
            _bonded = sorted(conect_map[_serial])
            _f.write(write_conect_line(_serial, _bonded))
        _f.write("END\n")

    logger.info(
        "Wrote assembled PDB: %s (chain A=protein, C=polymers, resid 1-%d)",
        assembled_path,
        chain_c_resid - 1,
    )
    return (assembled_path,)


@app.cell
def _(ForceField, PROTEIN_FF, SMALL_MOL_FF, conjugate_off, logger):
    """Cell 8: Create Interchange from the merged conjugate molecule."""

    conjugate_top = conjugate_off.to_topology()

    ff = ForceField(PROTEIN_FF, SMALL_MOL_FF)

    # Suppress per-atom INFO logging from OpenFF's LibraryCharges handler.
    import logging as _logging

    _nonbonded_logger = _logging.getLogger("openff.interchange.smirnoff._nonbonded")
    _prev_level = _nonbonded_logger.level
    _nonbonded_logger.setLevel(_logging.WARNING)
    try:
        interchange = ff.create_interchange(conjugate_top, charge_from_molecules=[conjugate_off])
    finally:
        _nonbonded_logger.setLevel(_prev_level)
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
    component_metadata,
    conjugated_sequences,
    conjugate_off,
    free_placement_results,
    free_polymer_offs,
    free_sequences,
    interchange,
    logger,
    placement_results: "list[PlacementResult]",
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
                    nbr = bond.atom1_index if bond.atom2_index == idx else bond.atom2_index
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

    mobile = _atoms_within_n_bonds(conjugate_off, linkage_seeds, LINKAGE_NEIGHBORHOOD_BONDS)
    logger.info("Mobile atoms (within %d bonds): %d", LINKAGE_NEIGHBORHOOD_BONDS, len(mobile))

    system = interchange.to_openmm_system()
    topology_omm = interchange.to_openmm_topology()
    _coords = conjugate_off.conformers[0].m_as(unit.angstrom)

    restraint = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter("k", 1000.0 * omm_unit.kilojoule_per_mole / omm_unit.nanometer**2)
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

    output_path = MINIMIZED_PDB
    with open(output_path, "w", encoding="utf-8") as _f:
        _atom_index = 0
        _updated_atoms = 0
        _preserved_atoms = 0
        n_simulated = len(minimized_coords)
        for _line in _old_lines:
            if _line.startswith(("ATOM", "HETATM")):
                if _atom_index < n_simulated:
                    _x, _y, _z = minimized_coords[_atom_index]
                    _f.write(f"{_line[:30]}{_x:8.3f}{_y:8.3f}{_z:8.3f}{_line[54:]}")
                    _updated_atoms += 1
                else:
                    _f.write(_line)
                    _preserved_atoms += 1
                _atom_index += 1
            else:
                _f.write(_line)

    if _updated_atoms != n_simulated:
        raise RuntimeError(
            f"PDB write mismatch: wrote coords for {_updated_atoms} atoms "
            f"but simulation had {n_simulated} atoms"
        )

    logger.info(
        "Minimized PDB coordinate update: %d simulated atoms updated, %d template atoms preserved",
        _updated_atoms,
        _preserved_atoms,
    )

    logger.info(
        "Minimization: E_before=%.2f, E_after=%.2f kJ/mol",
        energy_before,
        energy_after,
    )
    return (
        component_metadata,
        conjugated_sequences,
        energy_after,
        energy_before,
        free_polymer_offs,
        free_sequences,
        output_path,
    )


@app.cell
def _(
    EQUILIBRATED_PDB,
    PROTEIN_HEAVY_RESTRAINT_K,
    PROTEIN_PDB,
    PROTEIN_RMSD_THRESHOLD_A,
    VACUUM_EQ_FRICTION_PER_PS,
    VACUUM_EQ_STEPS,
    VACUUM_EQ_TEMP_K,
    VACUUM_EQ_TIMESTEP_FS,
    assembled_path,
    component_metadata,
    conjugate_off,
    interchange,
    logger,
    np,
    output_path,
    unit,
):
    """Cell 9.5: Protein-restrained vacuum equilibration + RMSD verification.

    After Cell 9's local linkage minimization, this cell performs:
    1. Protein-restrained energy minimization (all polymer atoms free)
    2. Short NVT equilibration at 310K (protein heavy atoms restrained)
    3. RMSD verification of protein vs original crystal structure

    The protein is held in place by strong harmonic restraints on all
    heavy atoms, while polymers (conjugated + free) relax freely.
    """
    try:
        import openmm
        from openmm import app as openmm_app
        from openmm import unit as omm_unit
    except ImportError as _exc:
        raise ImportError("OpenMM required for equilibration") from _exc

    # ── Platform selection: CUDA → OpenCL → CPU ───────────────────────────
    def _select_platform():
        """Select the fastest available OpenMM platform."""
        for name in ("CUDA", "OpenCL", "CPU"):
            try:
                platform = openmm.Platform.getPlatformByName(name)
                logger.info("Selected OpenMM platform: %s", name)
                return platform
            except Exception:
                continue
        raise RuntimeError("No suitable OpenMM platform found")

    eq_platform = _select_platform()

    # ── Build system from interchange ──────────────────────────────────────
    eq_system = interchange.to_openmm_system()
    eq_topology = interchange.to_openmm_topology()

    # Get current coordinates (post-Cell 9 minimization)
    eq_coords_angstrom = conjugate_off.conformers[0].m_as(unit.angstrom)
    eq_coords_nm = eq_coords_angstrom * 0.1

    # ── Add protein heavy-atom restraints ──────────────────────────────────
    protein_heavy_indices = component_metadata["restraint_groups"]["protein_heavy"]
    logger.info(
        "Applying %.0f kJ/mol/nm^2 restraints to %d protein heavy atoms",
        PROTEIN_HEAVY_RESTRAINT_K,
        len(protein_heavy_indices),
    )

    eq_restraint = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    eq_restraint.addGlobalParameter(
        "k", PROTEIN_HEAVY_RESTRAINT_K * omm_unit.kilojoule_per_mole / omm_unit.nanometer**2
    )
    eq_restraint.addPerParticleParameter("x0")
    eq_restraint.addPerParticleParameter("y0")
    eq_restraint.addPerParticleParameter("z0")

    for _idx in protein_heavy_indices:
        x0, y0, z0 = eq_coords_nm[_idx]
        eq_restraint.addParticle(_idx, [float(x0), float(y0), float(z0)])
    eq_system.addForce(eq_restraint)

    # ── Stage A: Protein-restrained minimization ───────────────────────────
    logger.info("Stage A: Protein-restrained energy minimization...")
    eq_integrator_min = openmm.VerletIntegrator(0.001 * omm_unit.picoseconds)
    eq_sim_min = openmm_app.Simulation(eq_topology, eq_system, eq_integrator_min, eq_platform)
    eq_sim_min.context.setPositions(eq_coords_nm * omm_unit.nanometers)

    state_before_min = eq_sim_min.context.getState(getEnergy=True)
    eq_energy_before_min = state_before_min.getPotentialEnergy().value_in_unit(
        omm_unit.kilojoule_per_mole
    )
    openmm.LocalEnergyMinimizer.minimize(eq_sim_min.context, tolerance=10.0, maxIterations=1000)
    state_after_min = eq_sim_min.context.getState(getEnergy=True, getPositions=True)
    eq_energy_after_min = state_after_min.getPotentialEnergy().value_in_unit(
        omm_unit.kilojoule_per_mole
    )
    minimized_positions = state_after_min.getPositions(asNumpy=True)
    logger.info(
        "  Minimization: E_before=%.2f, E_after=%.2f kJ/mol",
        eq_energy_before_min,
        eq_energy_after_min,
    )

    # ── Stage B: Short NVT equilibration ───────────────────────────────────
    logger.info(
        "Stage B: NVT equilibration at %.1f K for %d steps (%.1f ps)...",
        VACUUM_EQ_TEMP_K,
        VACUUM_EQ_STEPS,
        VACUUM_EQ_STEPS * VACUUM_EQ_TIMESTEP_FS / 1000.0,
    )

    # Need a fresh system since we can't change integrator on existing context
    eq_system_nvt = interchange.to_openmm_system()
    # Re-add the same restraint force
    eq_restraint_nvt = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    eq_restraint_nvt.addGlobalParameter(
        "k", PROTEIN_HEAVY_RESTRAINT_K * omm_unit.kilojoule_per_mole / omm_unit.nanometer**2
    )
    eq_restraint_nvt.addPerParticleParameter("x0")
    eq_restraint_nvt.addPerParticleParameter("y0")
    eq_restraint_nvt.addPerParticleParameter("z0")
    # Use the minimized positions as restraint reference
    min_pos_nm = minimized_positions.value_in_unit(omm_unit.nanometer)
    for _idx in protein_heavy_indices:
        x0, y0, z0 = min_pos_nm[_idx]
        eq_restraint_nvt.addParticle(_idx, [float(x0), float(y0), float(z0)])
    eq_system_nvt.addForce(eq_restraint_nvt)

    eq_integrator_nvt = openmm.LangevinMiddleIntegrator(
        VACUUM_EQ_TEMP_K * omm_unit.kelvin,
        VACUUM_EQ_FRICTION_PER_PS / omm_unit.picosecond,
        VACUUM_EQ_TIMESTEP_FS * omm_unit.femtosecond,
    )

    eq_sim_nvt = openmm_app.Simulation(eq_topology, eq_system_nvt, eq_integrator_nvt, eq_platform)
    eq_sim_nvt.context.setPositions(minimized_positions)
    eq_sim_nvt.context.setVelocitiesToTemperature(VACUUM_EQ_TEMP_K * omm_unit.kelvin)

    state_before_nvt = eq_sim_nvt.context.getState(getEnergy=True)
    eq_energy_before_nvt = state_before_nvt.getPotentialEnergy().value_in_unit(
        omm_unit.kilojoule_per_mole
    )

    eq_sim_nvt.step(VACUUM_EQ_STEPS)

    state_after_nvt = eq_sim_nvt.context.getState(getEnergy=True, getPositions=True)
    eq_energy_after_nvt = state_after_nvt.getPotentialEnergy().value_in_unit(
        omm_unit.kilojoule_per_mole
    )
    equilibrated_positions = state_after_nvt.getPositions(asNumpy=True)
    equilibrated_coords_angstrom = equilibrated_positions.value_in_unit(omm_unit.angstrom)

    logger.info(
        "  NVT: E_before=%.2f, E_after=%.2f kJ/mol",
        eq_energy_before_nvt,
        eq_energy_after_nvt,
    )

    # ── Write equilibrated PDB ─────────────────────────────────────────────
    with open(output_path, "r", encoding="utf-8") as _f:
        _template_lines = _f.readlines()

    equilibrated_output_path = EQUILIBRATED_PDB
    with open(equilibrated_output_path, "w", encoding="utf-8") as _f:
        _atom_index = 0
        _updated_atoms = 0
        _preserved_atoms = 0
        n_simulated = len(equilibrated_coords_angstrom)
        for _line in _template_lines:
            if _line.startswith(("ATOM", "HETATM")):
                if _atom_index < n_simulated:
                    _x, _y, _z = equilibrated_coords_angstrom[_atom_index]
                    _f.write(f"{_line[:30]}{_x:8.3f}{_y:8.3f}{_z:8.3f}{_line[54:]}")
                    _updated_atoms += 1
                else:
                    _f.write(_line)
                    _preserved_atoms += 1
                _atom_index += 1
            else:
                _f.write(_line)

    if _updated_atoms != n_simulated:
        raise RuntimeError(
            f"PDB write mismatch: wrote coords for {_updated_atoms} atoms "
            f"but simulation had {n_simulated} atoms"
        )

    logger.info("Wrote equilibrated PDB: %s", equilibrated_output_path)
    logger.info(
        "Equilibrated PDB coordinate update: %d simulated atoms updated, %d template atoms preserved",
        _updated_atoms,
        _preserved_atoms,
    )

    # ── RMSD verification ──────────────────────────────────────────────────
    # Compare protein heavy atoms in equilibrated structure vs original crystal.
    # Load original crystal structure coordinates.
    from rdkit import Chem as _Chem

    _crystal_rd = _Chem.MolFromPDBFile(PROTEIN_PDB, removeHs=False, sanitize=True)
    if _crystal_rd is None:
        raise ValueError("Failed to load crystal structure for RMSD comparison")
    _crystal_conf = _crystal_rd.GetConformer()
    _crystal_all_coords = np.array(
        [list(_crystal_conf.GetAtomPosition(i)) for i in range(_crystal_rd.GetNumAtoms())]
    )

    # protein_heavy_indices maps into the merged conjugate molecule.
    # The protein atoms in the merged molecule have the same order as the
    # original protein, minus the removed HZ atoms. We need to figure out
    # which original crystal atoms correspond to which merged heavy atoms.
    #
    # component_metadata["protein"]["heavy_atom_indices"] gives us the
    # indices in the merged molecule. These correspond to the non-removed
    # heavy atoms from the original protein, in order.
    #
    # For RMSD, we extract:
    #   - equilibrated_coords_angstrom[protein_heavy_indices] = equilibrated heavy
    #   - crystal heavy atom coords (skipping removed atoms, heavy only)
    #
    # Build the crystal heavy atom coords matching the merged ordering:
    # The merged protein was built by removing specific H atoms (hz_indices[1:])
    # from the original. So the protein atoms in merged order are all original
    # atoms except the removed ones, in original order. Heavy atoms are the
    # subset with atomic_number > 1.
    #
    # Since we don't have prot_removed_indices in this cell, we use a simpler
    # approach: the crystal heavy atom positions are at the same indices as
    # protein_heavy_indices but in the ORIGINAL protein (before any removal).
    # Actually, protein_heavy_indices are indices in the MERGED molecule.
    # We need the original indices.
    #
    # Simpler approach: extract backbone heavy (CA, C, N, O) from crystal by
    # PDB atom name, and compare to backbone heavy in the equilibrated structure.

    protein_backbone_heavy_indices = component_metadata["protein"]["backbone_heavy_atom_indices"]

    # Extract backbone coords from equilibrated structure
    eq_backbone_coords = equilibrated_coords_angstrom[protein_backbone_heavy_indices]

    # Extract backbone coords from crystal structure by matching PDB names
    _crystal_backbone_coords = []
    for _ai in range(_crystal_rd.GetNumAtoms()):
        _rd_atom = _crystal_rd.GetAtomWithIdx(_ai)
        _info = _rd_atom.GetPDBResidueInfo()
        if _info is None:
            continue
        _aname = _info.GetName().strip()
        if _aname in ("CA", "C", "N", "O") and _rd_atom.GetAtomicNum() > 1:
            _crystal_backbone_coords.append(_crystal_all_coords[_ai])
    _crystal_backbone_coords = np.array(_crystal_backbone_coords)

    crystal_backbone_count = len(_crystal_backbone_coords)
    equilibrated_backbone_count = len(eq_backbone_coords)
    if crystal_backbone_count != equilibrated_backbone_count:
        raise RuntimeError(
            "Backbone atom count mismatch for RMSD comparison: "
            f"crystal={crystal_backbone_count}, equilibrated={equilibrated_backbone_count}"
        )

    # Kabsch-aligned RMSD
    def _kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
        """Compute RMSD after optimal superposition (Kabsch algorithm)."""
        assert P.shape == Q.shape, f"Shape mismatch: {P.shape} vs {Q.shape}"
        centroid_p = P.mean(axis=0)
        centroid_q = Q.mean(axis=0)
        P_c = P - centroid_p
        Q_c = Q - centroid_q
        H = P_c.T @ Q_c
        U, _, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        S = np.diag([1.0, 1.0, d])
        R = Vt.T @ S @ U.T
        Q_aligned = (R @ P_c.T).T
        diff = Q_aligned - Q_c
        return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

    # Compute backbone RMSD
    protein_backbone_rmsd = _kabsch_rmsd(_crystal_backbone_coords, eq_backbone_coords)

    logger.info(
        "Protein backbone RMSD (vs crystal): %.4f A (threshold: %.1f A)",
        protein_backbone_rmsd,
        PROTEIN_RMSD_THRESHOLD_A,
    )

    if protein_backbone_rmsd > PROTEIN_RMSD_THRESHOLD_A:
        raise RuntimeError(
            "RMSD CHECK FAILED: protein backbone RMSD "
            f"{protein_backbone_rmsd:.4f} A exceeds threshold "
            f"{PROTEIN_RMSD_THRESHOLD_A:.1f} A"
        )

    protein_rmsd_passed = True
    logger.info("RMSD CHECK PASSED: protein structure preserved during equilibration.")

    return (
        eq_energy_after_min,
        eq_energy_after_nvt,
        eq_energy_before_min,
        eq_energy_before_nvt,
        equilibrated_output_path,
        protein_backbone_rmsd,
        protein_rmsd_passed,
    )


@app.cell
def _(
    energy_after,
    energy_before,
    eq_energy_after_min,
    eq_energy_after_nvt,
    eq_energy_before_min,
    eq_energy_before_nvt,
    equilibrated_output_path,
    logger,
    output_path,
    placement_results: "list[PlacementResult]",
    protein_backbone_rmsd,
    protein_rmsd_passed,
):
    """Cell 10: Summary and visualization."""
    import marimo as mo

    summary_rows = []
    for _i, _p in enumerate(placement_results):
        summary_rows.append(
            {
                "id": _i,
                "resid": _p.site.resid,
                "clearance_A": round(_p.min_protein_distance, 3),
            }
        )

    logger.info("Summary: %s", summary_rows)

    # Use equilibrated PDB for visualization if available
    _viz_path = equilibrated_output_path if equilibrated_output_path else output_path

    try:
        import py3Dmol

        with open(_viz_path, "r", encoding="utf-8") as _f:
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
        viz = mo.md(f"py3Dmol not available. PDB: `{_viz_path}`")

    summary = {
        "n_conjugates": len(placement_results),
        "linkage_minimization": {
            "energy_before_kj_mol": float(energy_before),
            "energy_after_kj_mol": float(energy_after),
        },
        "vacuum_equilibration": {
            "restraint_min_before_kj_mol": float(eq_energy_before_min),
            "restraint_min_after_kj_mol": float(eq_energy_after_min),
            "nvt_before_kj_mol": float(eq_energy_before_nvt),
            "nvt_after_kj_mol": float(eq_energy_after_nvt),
        },
        "protein_backbone_rmsd_A": float(protein_backbone_rmsd),
        "protein_rmsd_passed": protein_rmsd_passed,
        "output_pdb": str(equilibrated_output_path),
        "rows": summary_rows,
    }
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
