"""Runtime validation helpers for simulation configuration references."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def collect_reference_warnings(config: Any) -> list[str]:
    """Return warnings for missing files referenced by a simulation config.

    These checks intentionally live outside the Pydantic schema so lightweight
    configuration parsing can succeed before external structure files are staged.

    Parameters
    ----------
    config : Any
        Simulation configuration or config-like object with PolyzyMD schema
        attributes.

    Returns
    -------
    list[str]
        Human-readable warnings for missing referenced files or directories.
    """

    warnings: list[str] = []
    enzyme = getattr(config, "enzyme", None)
    _check_file(warnings, getattr(enzyme, "pdb_path", None), "enzyme PDB")

    substrate = getattr(config, "substrate", None)
    if substrate is not None:
        _check_file(warnings, getattr(substrate, "sdf_path", None), "substrate SDF")

    polymers = getattr(config, "polymers", None)
    if polymers is not None and bool(getattr(polymers, "enabled", False)):
        _check_polymer_references(warnings, polymers)

    return warnings


def _check_file(warnings: list[str], path_value: Any, label: str) -> None:
    """Append a warning when a referenced file path is missing."""

    if path_value is None:
        return
    path = Path(path_value)
    if not path.is_file():
        warnings.append(f"Missing {label}: {path}")


def _check_directory(warnings: list[str], path_value: Any, label: str) -> Path | None:
    """Append a warning when a referenced directory path is missing."""

    if path_value is None:
        return None
    path = Path(path_value)
    if not path.is_dir():
        warnings.append(f"Missing {label}: {path}")
        return None
    return path


def _check_polymer_references(warnings: list[str], polymers: Any) -> None:
    """Check polymer SDF and reaction-template references."""

    generation_mode = str(getattr(polymers, "generation_mode", "")).lower()
    if generation_mode.endswith("cached"):
        sdf_directory = _check_directory(
            warnings,
            getattr(polymers, "sdf_directory", None),
            "polymer SDF directory",
        )
        if sdf_directory is not None:
            _check_cached_polymer_sdfs(warnings, polymers, sdf_directory)

    reactions = getattr(polymers, "reactions", None)
    if reactions is not None:
        for field_name in ("initiation", "polymerization", "termination"):
            _check_file(
                warnings,
                getattr(reactions, field_name, None),
                f"polymer {field_name} reaction template",
            )


def _check_cached_polymer_sdfs(warnings: list[str], polymers: Any, sdf_directory: Path) -> None:
    """Warn when a cached-polymer directory has no matching SDF files."""

    type_prefix = getattr(polymers, "type_prefix", None)
    length = getattr(polymers, "length", None)
    if not type_prefix or length is None:
        return

    pattern = f"{type_prefix}_seq=*_{length}-mer_charged.sdf"
    if not any(sdf_directory.glob(pattern)):
        warnings.append(
            "Missing cached polymer SDF files: " f"no files matching {sdf_directory / pattern}"
        )


#: Elements the NAGL charge model handles (the atom features of
#: ``openff-gnn-am1bcc-0.1.0-rc.3.pt``).
NAGL_ELEMENTS = frozenset({"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"})


def require_inputs(config: Any) -> None:
    """Raise ``ValueError`` when an input of a new build is missing or cannot be used.

    Cheap checks of what the build reads: the enzyme PDB has atoms, the
    substrate SDF holds the configured conformer, co-solvent and monomer SMILES
    parse, monomers carry the group the initiation reaction needs, the charge
    method can charge each molecule, and the force fields are installed.
    :meth:`SimulationConfig.require_buildable` calls this.
    """
    pdb = Path(config.enzyme.pdb_path)
    if not pdb.is_file():
        raise ValueError(f"enzyme PDB {pdb} does not exist.\nfix: correct enzyme.pdb_path.")
    with open(pdb, errors="replace") as stream:
        if not any(line.startswith(("ATOM", "HETATM")) for line in stream):
            raise ValueError(
                f"enzyme PDB {pdb} has no ATOM or HETATM records.\n"
                "fix: give a PDB file of the prepared protein (polyzymd clean-pdb writes one)."
            )

    from rdkit import Chem

    substrate = config.substrate
    if substrate is not None:
        sdf = Path(substrate.sdf_path)
        if not sdf.is_file():
            raise ValueError(
                f"substrate SDF {sdf} does not exist.\nfix: correct substrate.sdf_path."
            )
        molecules = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
        if substrate.conformer_index >= len(molecules):
            raise ValueError(
                f"substrate SDF {sdf} holds {len(molecules)} conformer(s), so conformer_index "
                f"{substrate.conformer_index} does not exist.\n"
                f"fix: set conformer_index between 0 and {max(len(molecules) - 1, 0)}."
            )
        molecule = molecules[substrate.conformer_index]
        if molecule is None:
            raise ValueError(
                f"RDKit cannot read conformer {substrate.conformer_index} of {sdf}.\n"
                "fix: write the SDF again from your docking or drawing program."
            )
        _require_chargeable(molecule, substrate.charge_method, f"substrate {substrate.name!r}")

    from polyzymd.data.solvent_molecules import is_bundled_solvent, split_counter_ions

    for cosolvent in config.solvent.co_solvents:
        if is_bundled_solvent(cosolvent.name):
            continue
        what = f"co-solvent {cosolvent.name!r}"
        smiles, _, _ = split_counter_ions(cosolvent.smiles, cosolvent.name)
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(
                f"{what}: RDKit cannot read the SMILES {cosolvent.smiles!r}.\n"
                "fix: correct the SMILES."
            )
        if molecule.GetNumAtoms() == 1 and molecule.GetAtomWithIdx(0).GetFormalCharge():
            raise ValueError(
                f"{what} is the ion {smiles}; the build adds only Na+ and Cl- ions.\n"
                "fix: remove this co-solvent and use solvent.ions.nacl_concentration."
            )
        _require_chargeable(molecule, cosolvent.charge_method, what)

    polymers = config.polymers
    if polymers is not None and polymers.enabled and polymers.generation_mode.value == "dynamic":
        from rdkit.Chem import AllChem

        initiation = Path(polymers.reactions.initiation)
        if not initiation.is_file():
            raise ValueError(
                f"polymer initiation reaction {initiation} does not exist.\n"
                "fix: correct polymers.reactions.initiation, or set it to default."
            )
        reaction = AllChem.ReactionFromRxnFile(str(initiation))
        templates = [
            reaction.GetReactantTemplate(i) for i in range(reaction.GetNumReactantTemplates())
        ]
        for monomer in polymers.monomers:
            what = f"monomer {monomer.label!r}"
            molecule = Chem.MolFromSmiles(monomer.smiles)
            if molecule is None:
                raise ValueError(
                    f"{what}: RDKit cannot read the SMILES {monomer.smiles!r}.\n"
                    "fix: correct the SMILES."
                )
            molecule = Chem.AddHs(molecule)
            if not any(molecule.HasSubstructMatch(template) for template in templates):
                raise ValueError(
                    f"{what} ({monomer.smiles}) has no group the initiation reaction {initiation.name} "
                    "reacts with, so no polymer can be grown from it.\n"
                    "fix: give the monomer with its polymerisable group (a methacrylate for the "
                    "default reactions), or reactions that match it."
                )
            _require_chargeable(molecule, polymers.charger, what)

    from openff.toolkit.typing.engines.smirnoff.forcefield import get_available_force_fields

    available = set(get_available_force_fields())
    for key in ("protein", "small_molecule"):
        name = getattr(config.force_field, key)
        if name not in available and not Path(name).is_file():
            raise ValueError(
                f"force_field.{key} {name!r} is not an installed force field.\n"
                "fix: correct the name (the default is "
                f"{type(config.force_field).model_fields[key].default}), or give the path of an "
                ".offxml file."
            )


def _require_chargeable(molecule: Any, method: Any, what: str) -> None:
    """Raise ``ValueError`` when ``method`` cannot assign charges to an RDKit molecule."""
    method = getattr(method, "value", method)
    if method == "nagl":
        elements = {atom.GetSymbol() for atom in molecule.GetAtoms()}
        others = sorted(elements - NAGL_ELEMENTS)
        if others:
            raise ValueError(
                f"{what} contains {', '.join(others)}, which the NAGL charge model does not "
                f"cover (it covers {', '.join(sorted(NAGL_ELEMENTS))}).\n"
                "fix: set charge_method: am1bcc (needs AmberTools), or leave the molecule out."
            )
    elif method == "am1bcc":
        from openff.toolkit.utils import AmberToolsToolkitWrapper, OpenEyeToolkitWrapper

        if not (AmberToolsToolkitWrapper.is_available() or OpenEyeToolkitWrapper.is_available()):
            raise ValueError(
                f"{what}: charge_method am1bcc needs AmberTools or OpenEye, and neither is "
                "installed in this environment.\nfix: use charge_method: nagl."
            )
    elif method == "espaloma":
        import importlib.util

        if importlib.util.find_spec("espaloma_charge") is None:
            raise ValueError(
                f"{what}: charge_method espaloma needs the espaloma-charge package, which is not "
                "installed in this environment.\nfix: use charge_method: nagl."
            )
