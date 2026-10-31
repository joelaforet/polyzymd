"""Runtime validation helpers for simulation configuration references."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence


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

    _check_am1bcc(warnings, config)
    return warnings


def stage_temperature_mismatch(stage: Any, temperature: float, *, last: bool) -> float | None:
    """Return the stage temperature that differs from the production temperature.

    A constant-temperature stage is compared by its ``temperature``. A ramp is
    compared by its ``temperature_end``, and only when it is the last stage,
    because production starts where that ramp ends.

    Parameters
    ----------
    stage : Any
        An equilibration stage.
    temperature : float
        The production temperature, ``thermodynamics.temperature`` (K).
    last : bool
        Whether the stage is the last one before production.

    Returns
    -------
    float | None
        The stage temperature (K) when it differs, else ``None``.
    """
    if getattr(stage, "is_temperature_ramping", False):
        value = stage.temperature_end if last else None
    else:
        value = getattr(stage, "temperature", None)
    if value is None or abs(float(value) - float(temperature)) < 1e-6:
        return None
    return float(value)


def _stage_start(stage: Any) -> float | None:
    """The temperature (K) a stage starts at: a ramp's start, else its temperature."""
    if getattr(stage, "is_temperature_ramping", False):
        return stage.temperature_start
    return getattr(stage, "temperature", None)


def stage_temperature_warning(stages: Sequence[Any], index: int, temperature: float) -> str | None:
    """Return the warning for one equilibration stage's temperature, or ``None``.

    A constant-temperature stage, and a ramp that is the last stage, are
    compared with ``temperature`` (production's). A ramp before another stage
    is compared with the temperature that stage starts at, since a mismatch is
    a sudden jump between the two.

    Parameters
    ----------
    stages : Sequence[Any]
        The equilibration stages, in order.
    index : int
        The stage to check.
    temperature : float
        The production temperature, ``thermodynamics.temperature`` (K).

    Returns
    -------
    str | None
        The warning, or ``None`` when the temperature matches.
    """
    stage = stages[index]
    last = index == len(stages) - 1
    ramp = getattr(stage, "is_temperature_ramping", False)
    if ramp and not last:
        following = stages[index + 1]
        start = _stage_start(following)
        end = stage.temperature_end
        if start is None or end is None or abs(float(end) - float(start)) < 1e-6:
            return None
        return (
            f"equilibration stage '{stage.name}' ends at {float(end):g} K, but the next stage "
            f"'{following.name}' starts at {float(start):g} K, a sudden jump of "
            f"{abs(float(start) - float(end)):g} K. If you changed the temperature of this "
            f"condition, also set simulation_phases.equilibration_stages.{index}.temperature_end."
        )
    value = stage_temperature_mismatch(stage, temperature, last=last)
    if value is None:
        return None
    key = "temperature_end" if ramp else "temperature"
    where = "ends at" if ramp else "runs at"
    return (
        f"equilibration stage '{stage.name}' {where} {value:g} K, but production runs at "
        f"thermodynamics.temperature = {temperature:g} K. If you changed the temperature "
        f"of this condition, also set simulation_phases.equilibration_stages.{index}.{key}."
    )


def collect_protocol_warnings(config: Any) -> list[str]:
    """Return warnings for equilibration temperatures that do not join up, and for
    restraint selections of a modified protein that pick residues by number.

    Production runs at ``thermodynamics.temperature``. An equilibration stage
    at another temperature is allowed, so this is a warning, not an error: a
    condition copied to a new temperature often keeps the old stage
    temperatures by mistake. A heating ramp is checked against the stage after
    it, or against production when it is the last stage.

    Parameters
    ----------
    config : Any
        Simulation configuration.

    Returns
    -------
    list[str]
        One warning for each stage whose temperature does not join up.
    """
    warnings = _modification_selection_warnings(config)
    phases = getattr(config, "simulation_phases", None)
    thermodynamics = getattr(config, "thermodynamics", None)
    if thermodynamics is None or phases is None:
        return warnings
    temperature = float(thermodynamics.temperature)
    stages = list(getattr(phases, "equilibration_stages", None) or [])
    found = (stage_temperature_warning(stages, i, temperature) for i in range(len(stages)))
    return [warning for warning in found if warning is not None] + warnings


def _modification_selection_warnings(config: Any) -> list[str]:
    """Restraint atoms on a protein with attachments, selected by residue number.

    The selection is matched among all the protein component's atoms, its
    attachment segments included, whose residues are numbered from 1.
    """
    from polyzymd.config.modifications import has_attachments, unscoped_residue_selection

    modified = {
        s.name for s in getattr(config, "solutes", None) or [] if has_attachments(s)
    }
    warnings = []
    for restraint in getattr(config, "restraints", None) or []:
        for which in ("atom1", "atom2"):
            atom = getattr(restraint, which)
            if atom.component in modified and unscoped_residue_selection(atom.selection):
                warnings.append(
                    f"restraint {restraint.name} {which} selects '{atom.selection}' in "
                    f"{atom.component}, whose attachments (glycans, grafts, labels) are "
                    "segments numbered from resid 1, so the selection may also match them; "
                    "the build's restraint.atom_identity check names the atom it picked."
                )
    return warnings


def _check_am1bcc(warnings: list[str], config: Any) -> None:
    """Warn when the config asks for AM1-BCC charges and no toolkit here computes them."""
    from polyzymd.utils import charging

    substrate = getattr(config, "substrate", None)
    solvent = getattr(config, "solvent", None)
    polymers = getattr(config, "polymers", None)
    methods = [getattr(substrate, "charge_method", None)]
    methods += [cs.charges.method for cs in getattr(solvent, "co_solvents", None) or []]
    if polymers is not None and getattr(polymers, "enabled", False):
        methods.append(getattr(polymers, "charger", None))
    if "am1bcc" in {getattr(m, "value", m) for m in methods} and not charging.am1bcc_available():
        warnings.append(charging.AM1BCC_NEEDS_AMBERTOOLS)


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
    """Check the SDF directory of cached polymers."""

    generation_mode = str(getattr(polymers, "generation_mode", "")).lower()
    if generation_mode.endswith("cached"):
        sdf_directory = _check_directory(
            warnings,
            getattr(polymers, "sdf_directory", None),
            "polymer SDF directory",
        )
        if sdf_directory is not None:
            _check_cached_polymer_sdfs(warnings, polymers, sdf_directory)


def _check_cached_polymer_sdfs(warnings: list[str], polymers: Any, sdf_directory: Path) -> None:
    """Warn when a cached-polymer directory has no matching SDF files."""

    type_prefix = getattr(polymers, "type_prefix", None)
    length = getattr(polymers, "length", None)
    if not type_prefix or length is None:
        return

    pattern = f"{type_prefix}_seq=*_{length}-mer_charged.sdf"
    if not any(sdf_directory.glob(pattern)):
        warnings.append(
            f"Missing cached polymer SDF files: no files matching {sdf_directory / pattern}"
        )


#: Elements the NAGL charge model handles (the atom features of
#: ``openff-gnn-am1bcc-0.1.0-rc.3.pt``).
NAGL_ELEMENTS = frozenset({"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"})


def _key(config: Any, solute: Any, field: str) -> str:
    """The config key of ``field`` on ``solute``, as the loader names it (``solutes.0.pdb_path``)."""
    index = next(i for i, entry in enumerate(config.solutes) if entry is solute)
    return f"solutes.{index}.{field}"


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
        raise ValueError(
            f"protein PDB {pdb} does not exist.\n"
            f"fix: correct {_key(config, config.enzyme, 'pdb_path')}."
        )
    with open(pdb, errors="replace") as stream:
        if not any(line.startswith(("ATOM", "HETATM")) for line in stream):
            raise ValueError(
                f"protein PDB {pdb} has no ATOM or HETATM records.\n"
                "fix: give a PDB file of the prepared protein (polyzymd clean-pdb writes one)."
            )

    from rdkit import Chem

    substrate = config.substrate
    if substrate is not None:
        sdf = Path(substrate.sdf_path)
        if not sdf.is_file():
            raise ValueError(
                f"small-molecule SDF {sdf} does not exist.\n"
                f"fix: correct {_key(config, substrate, 'sdf_path')}."
            )
        molecules = list(Chem.SDMolSupplier(str(sdf), removeHs=False))
        if substrate.conformer_index >= len(molecules):
            raise ValueError(
                f"small-molecule SDF {sdf} holds {len(molecules)} conformer(s), so conformer_index "
                f"{substrate.conformer_index} does not exist.\n"
                f"fix: set conformer_index between 0 and {max(len(molecules) - 1, 0)}."
            )
        molecule = molecules[substrate.conformer_index]
        if molecule is None:
            raise ValueError(
                f"RDKit cannot read conformer {substrate.conformer_index} of {sdf}.\n"
                "fix: write the SDF again from your docking or drawing program."
            )
        _require_chargeable(molecule, substrate.charge_method, f"solute {substrate.name!r}")
        _require_restraint_atom_names(config, substrate, molecule)

    from polyzymd.data.solvent_molecules import is_bundled_solvent, split_counter_ions

    for cosolvent in config.solvent.co_solvents:
        if cosolvent.charges.method == "preset" and is_bundled_solvent(cosolvent.name):
            _require_shipped_molecule(cosolvent)
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
        _require_chargeable(molecule, cosolvent.charges.method, what)

    polymers = config.polymers
    if polymers is not None and polymers.enabled and polymers.generation_mode.value == "dynamic":
        from polyzymd.builders.polymer_generator import parse_star_smiles, require_mbuild_polymer

        units = [
            (f"monomer {m.label!r}", parse_star_smiles(m.smiles, 2, f"monomer {m.label!r}"))
            for m in polymers.monomers
        ]
        caps = [
            (f"end_groups.{end} {option.label!r}", option)
            for end in ("head", "tail")
            for option in getattr(polymers.end_groups, end)
        ]
        units += [(what, parse_star_smiles(option.smiles, 1, what)) for what, option in caps]
        for what, unit in units:
            _require_chargeable(unit.mol, polymers.charger, what)
        from polyzymd.config.schema import reserved_residue_names

        named = [(f"monomer {m.label!r}", m.residue_name) for m in polymers.monomers]
        named += [(what, option.residue_name) for what, option in caps]
        for what, residue_name in named:
            if residue_name and residue_name.upper() in reserved_residue_names():
                raise ValueError(
                    f"{what} has residue_name '{residue_name}', an amino-acid, nucleic-acid, "
                    "water or ion name, so selections such as 'protein' would include the "
                    "polymer.\n"
                    f"fix: set the residue_name of {what} in {_key(config, polymers, 'monomers')} "
                    "or end_groups to another 3-letter name."
                )
        try:
            require_mbuild_polymer()
        except ImportError as exc:
            raise ValueError(str(exc)) from exc

    if polymers is not None and polymers.enabled:
        _require_cached_chemistry(config, polymers)

    from openff.toolkit.typing.engines.smirnoff.forcefield import get_available_force_fields

    available = set(get_available_force_fields())
    named = [(f"solute {s.name!r}", s.force_field) for s in config.solutes]
    named.append(("solvent", config.solvent.force_field))
    for what, files in named:
        for name in files:
            if name not in available and not Path(name).is_file():
                raise ValueError(
                    f"{what}: force_field {name!r} is not an installed force field.\n"
                    "fix: correct the name (for example ff14sb_off_impropers_0.0.4.offxml or "
                    "openff-2.3.0.offxml), or give the path of an .offxml file."
                )

    for solute in config.solutes:
        if getattr(solute, "is_modified", False):
            _require_modification_inputs(config, solute, available)

    from polyzymd.utils.charging import installed_nagl_model

    nagl_model = getattr(config, "nagl_model", None)
    if nagl_model is not None and installed_nagl_model(nagl_model) is None:
        from openff.nagl_models import list_available_nagl_models

        installed = ", ".join(sorted(path.name for path in list_available_nagl_models()))
        raise ValueError(
            f"charges model {nagl_model!r} is not an installed NAGL model.\n"
            f"fix: in each charges: {{method: nagl, model: ...}}, name one of {installed}, "
            "or give the path of a model file."
        )


def _require_modification_inputs(config: Any, solute: Any, available: set) -> None:
    """Check the files, force fields and NAGL models a protein's modifications name.

    Every file must exist (CCD definitions, fragment PDBs, sites files, the
    glycan file of each accession), a sites CSV must have a ``site`` column
    (and ``glycan`` where it names glycans), every force field and NAGL model
    of ``modification_parameters`` must be installed, and no atom may be
    modified twice.
    """
    from polyzymd.config.modifications import GLYCAM_FORCE_FIELD, glycan_file

    def need(path: Any, what: str, key: str) -> Path:
        path = Path(path)
        if not path.is_file():
            raise ValueError(
                f"{what} {path} does not exist.\nfix: correct {_key(config, solute, key)}."
            )
        return path

    for n, mutation in enumerate(solute.mutations):
        ccd_file = getattr(mutation.to, "ccd_file", None)
        if ccd_file is not None:
            need(ccd_file, "CCD file", f"mutations.{n}.to.ccd_file")
        if mutation.sites_file is not None:
            _sites_csv(need(mutation.sites_file, "sites file", f"mutations.{n}.sites_file"), [])
    for n, entry in enumerate(solute.modifications):
        verb, arguments = entry.verb, entry.arguments
        where = f"modifications.{n}.{verb}"
        fragment = getattr(arguments, "fragment", None)
        if fragment is not None and fragment.pdb_path is not None:
            need(fragment.pdb_path, "fragment PDB", f"{where}.fragment.pdb_path")
        if fragment is not None and fragment.ccd_file is not None:
            need(fragment.ccd_file, "fragment CCD file", f"{where}.fragment.ccd_file")
        rows = []
        if arguments.sites_file is not None:
            wanted = ["glycan"] if verb == "attach" and fragment is None else []
            rows = _sites_csv(
                need(arguments.sites_file, "sites file", f"{where}.sites_file"), wanted
            )
        accessions = list(getattr(arguments, "glycans", None) or {}) + [
            row["glycan"] for row in rows if row.get("glycan")
        ]
        for accession in accessions:
            path = glycan_file(arguments.glycan_library, accession)
            if not path.is_file():
                raise ValueError(
                    f"no file for glycan {accession}: expected {path} (or "
                    f"{path.with_name(accession + '.pdb')}), with GLYCAM06 residue names and "
                    "CONECT records.\n"
                    f"fix: put the file in {path.parent}, or set {_key(config, solute, where)}"
                    ".glycan_library to the folder that holds it."
                )
    if solute.modification_parameters is not None:
        from polyzymd.utils.charging import installed_nagl_model

        for site_type, chosen in solute.modification_parameters.as_dict().items():
            key = _key(config, solute, f"modification_parameters.{site_type}")
            for name in chosen.force_field:
                if Path(name).name == GLYCAM_FORCE_FIELD or name in available:
                    continue
                if not Path(name).is_file():
                    raise ValueError(
                        f"{key}: force_field {name!r} is not an installed force field.\n"
                        "fix: name an installed OpenFF force field (for example "
                        f"openff-2.3.0.offxml), {GLYCAM_FORCE_FIELD} for glycans, or the path "
                        "of an .offxml file."
                    )
            model = chosen.charges.model
            if model is not None and installed_nagl_model(model) is None:
                raise ValueError(
                    f"{key}: charges model {model!r} is not an installed NAGL model.\n"
                    "fix: name an installed model (openff-gnn-am1bcc-1.0.0.pt) or the path of "
                    "a model file."
                )
    _require_one_change_per_atom(config, solute)


def _sites_csv(path: Path, wanted: list[str]) -> list[dict[str, str]]:
    """The rows of a sites CSV, refused when it lacks ``site`` or a ``wanted`` column."""
    import csv

    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = [name.strip() for name in reader.fieldnames or []]
        rows = [
            {key.strip(): (value or "").strip() for key, value in row.items() if key}
            for row in reader
        ]
    missing = [name for name in ["site", *wanted] if name not in columns]
    if missing:
        raise ValueError(
            f"sites file {path} has no {', '.join(repr(m) for m in missing)} column "
            f"(columns: {', '.join(columns) or 'none'}).\n"
            "fix: give a header row with `site` (A:N58) and, for glycans, `glycan` (an accession)."
        )
    if not rows:
        raise ValueError(f"sites file {path} has a header but no rows.\nfix: list the sites.")
    return rows


def _require_one_change_per_atom(config: Any, solute: Any) -> None:
    """Refuse two attachments at one atom, two identical protonation changes, or two
    mutations of one residue; the message names both entries."""
    from polyzymd.config.modifications import declared_modifications

    seen: dict[tuple, str] = {}
    for item in declared_modifications(
        solute.mutations, solute.modifications, solute.modification_parameters
    ):
        if item["kind"] == "mutate":
            key = ("mutate", item["site"])
        elif item["kind"] == "attach":
            key = ("attach", item["site"], item["atom_name"])
        else:
            key = (item["kind"], item["site"], item["atom_name"])
        if key in seen:
            what = " ".join(str(part) for part in key[1:])
            raise ValueError(
                f"{_key(config, solute, item['entry'])} and "
                f"{_key(config, solute, seen[key])} both {key[0]} {what}.\n"
                "fix: keep one of them; each residue is mutated once and each atom takes one "
                "attachment."
            )
        seen[key] = item["entry"]


def _require_cached_chemistry(config: Any, polymers: Any) -> None:
    """Raise ``ValueError`` when cached chains were built from other SMILES or caps.

    Only chains that dynamic generation made carry the build metadata this
    compares; the build makes the same check for the chains it loads.
    """
    if polymers.generation_mode.value != "cached" or polymers.sdf_directory is None:
        return
    if not any(m.smiles for m in polymers.monomers) and polymers.end_groups is None:
        return
    from polyzymd.builders.polymer import PolymerBuilder

    end_groups = None
    if polymers.end_groups is not None:
        end_groups = {
            end: [option.model_dump() for option in getattr(polymers.end_groups, end)]
            for end in ("head", "tail")
        }
    builder = PolymerBuilder(
        characters=[m.label for m in polymers.monomers],
        probabilities=[m.probability for m in polymers.monomers],
        length=polymers.length,
        type_prefix=polymers.type_prefix,
        sdf_directory=polymers.sdf_directory,
        monomer_smiles={(m.name or m.label): m.smiles for m in polymers.monomers if m.smiles},
        monomer_names={m.label: (m.name or m.label) for m in polymers.monomers},
        end_groups=end_groups,
    )
    directory = Path(polymers.sdf_directory)
    for sdf in sorted(directory.glob(f"{polymers.type_prefix}_seq=*_charged.sdf")):
        builder._check_cached_chemistry(sdf)


def small_molecule_atom_names(molecule: Any) -> list[str]:
    """Return the atom names a small molecule read from an SDF gets in the built system.

    OpenFF keeps the atom names of a molecule only when every atom has one and
    no two are alike. An SDF has no atom names unless it carries an
    ``atom.prop._Name`` data field, so OpenFF names each atom when it writes
    the topology: element symbol, count of that element so far in the file's
    atom order, and an ``x`` suffix (``C1x``, ``C2x``, ``O1x``, ``H1x``; see
    ``openff.toolkit.Molecule.generate_unique_atom_names``).

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        The molecule, with its hydrogens, in the SDF's atom order.
    """
    atoms = list(molecule.GetAtoms())
    given = [atom.GetProp("_Name") if atom.HasProp("_Name") else "" for atom in atoms]
    if len(set(given)) == len(given):
        return given
    counts: dict[str, int] = {}
    names = []
    for atom in atoms:
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
        names.append(f"{symbol}{counts[symbol]}x")
    return names


#: Words that end the list of values after ``name`` in a selection.
_SELECTION_STOP_WORDS = frozenset({"and", "or", "not", "to"})


def selection_atom_names(selection: str) -> list[str] | None:
    """Return the atom names a selection asks for with ``name``, or None when they cannot be listed.

    ``"resname LIG and name C1 C2"`` gives ``["C1", "C2"]``; a selection
    without ``name`` gives ``[]``. A regular expression (``name =~ 'C.*'``)
    or a ``not`` gives None: a missing name then need not mean that the
    selection matches nothing.
    """
    import re
    import shlex

    try:
        words = shlex.split(re.sub(r"([()])", r" \1 ", selection))
    except ValueError:
        return None
    if any(word.lower() == "not" for word in words):
        return None
    names: list[str] = []
    for position, word in enumerate(words):
        if word != "name":
            continue
        for value in words[position + 1 :]:
            if value == "=~":
                return None
            if value == "==":
                continue
            if value in ("(", ")") or value.lower() in _SELECTION_STOP_WORDS:
                break
            names.append(value)
    return names


def _require_restraint_atom_names(config: Any, solute: Any, molecule: Any) -> None:
    """Raise ``ValueError`` when a restraint names atoms the small molecule will not have.

    The built system names the atoms of a small molecule from an SDF without
    atom names ``C1x``, ``C2x``, ... (:func:`small_molecule_atom_names`), so
    ``name C1`` matches nothing and the build would fail only after solvation
    and parameterisation. Only selections whose every ``name`` value is
    missing are refused; any other selection is left to the build.
    """
    available = small_molecule_atom_names(molecule)
    if not all(available):
        return  # an unnamed single atom keeps its empty name; no name selects it
    for index, restraint in enumerate(getattr(config, "restraints", None) or []):
        if not restraint.enabled:
            continue
        for key in ("atom1", "atom2"):
            atom = getattr(restraint, key)
            if atom.component != solute.name:
                continue
            asked = selection_atom_names(atom.selection)
            if not asked or any(name in available for name in asked):
                continue
            raise ValueError(
                f"restraint {restraint.name!r}: {key} selection {atom.selection!r} names "
                f"{', '.join(asked)}, which the small molecule {solute.name!r} will not have. "
                f"In the built system its atoms, in the order of {Path(solute.sdf_path).name}, "
                f"are named {', '.join(available)}.\n"
                f"fix: in restraints.{index}.{key}.selection, use one of these names "
                f"(for example name {available[0]})."
            )


def _require_shipped_molecule(cosolvent: Any) -> None:
    """Raise ``ValueError`` when a preset co-solvent's SMILES is not the molecule PolyzyMD ships."""
    from rdkit import Chem

    from polyzymd.data.solvent_molecules import _SOLVENTS_DIR, split_counter_ions

    library = _SOLVENTS_DIR / f"{cosolvent.name.lower().strip()}.sdf"
    shipped = next(iter(Chem.SDMolSupplier(str(library), removeHs=False)), None)
    smiles, _, _ = split_counter_ions(cosolvent.smiles, cosolvent.name)
    given = Chem.MolFromSmiles(smiles)
    if shipped is None or given is None:
        return
    if Chem.MolToSmiles(Chem.RemoveHs(shipped)) != Chem.MolToSmiles(given):
        raise ValueError(
            f"co-solvent {cosolvent.name!r}: its SMILES {cosolvent.smiles!r} is not the "
            f"{cosolvent.name} PolyzyMD ships charges for, so charges preset cannot apply.\n"
            "fix: remove smiles to use the shipped molecule, or set charges: {method: nagl} "
            "(or am1bcc, espaloma)."
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
                "fix: set charges: {method: am1bcc} (needs AmberTools), or leave the molecule "
                "out."
            )
    elif method == "am1bcc":
        from polyzymd.utils import charging

        if not charging.am1bcc_available():
            raise ValueError(
                f"{what}: charges method am1bcc needs AmberTools or OpenEye, and neither is "
                "installed in this environment.\nfix: use charges: {method: nagl}."
            )
    elif method == "espaloma":
        import importlib.util

        if importlib.util.find_spec("espaloma_charge") is None:
            raise ValueError(
                f"{what}: charges method espaloma needs the espaloma-charge package, which is not "
                "installed in this environment.\nfix: use charges: {method: nagl}."
            )
