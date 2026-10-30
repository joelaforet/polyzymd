"""Build a protein solute that has ``mutations:`` or ``modifications:``.

The protein is loaded and changed through :mod:`polyzymd.modifying_proteins`
(mutations first, then the modifications in order, then the relaxation the
config asks for), written to the build folder, and parameterized by
``parameterize``: ff14SB outside the sites, each site from its type's
parameter set. The builder places the result like any protein; its terms and
charges replace the ones the combined force field gave the protein, GLYCAM06j
1-4 pairs are written into both engines' files, and its component record
carries one domain per site for the physics audit.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from openff.interchange import Interchange
    from openff.toolkit import Molecule

    from polyzymd.config.schema import ProteinSolute
    from polyzymd.core.components import ComponentRecord

LOGGER = logging.getLogger(__name__)

#: File the modified protein is written to, in the run folder.
MODIFIED_PROTEIN_PDB = "modified_protein.pdb"


def config_spec(
    solute: "ProteinSolute", replicate: int | None = None
) -> dict[str, list[dict[str, Any]]]:
    """Return the solute's blocks as the plain data ``apply_spec`` takes.

    A grafted chain without ``random_seed`` gets ``replicate``, which seeds its
    draw: each replicate then carries its own chain.
    """

    def plain(model) -> dict[str, Any]:
        return model.model_dump(mode="json", exclude_none=True, exclude_defaults=True)

    def modification(entry) -> dict[str, Any]:
        arguments = plain(entry.arguments)
        fragment = getattr(entry.arguments, "fragment", None)
        if fragment is not None and fragment.polymer is not None:
            chosen = solute.modification_parameters and solute.modification_parameters.polymer
            arguments["fragment"] = fragment.as_source(list(chosen.force_field) if chosen else [])
            chain = arguments["fragment"]["polymer"]
            if chain.get("random_seed") is None and replicate is not None:
                chain["replicate"] = replicate
        return {entry.verb: arguments}

    return {
        "mutations": [plain(mutation) for mutation in solute.mutations],
        "modifications": [modification(entry) for entry in solute.modifications],
    }


def parameter_sets(solute: "ProteinSolute") -> dict[str, Any]:
    """Return the solute's ``modification_parameters`` as ParameterSets, by site type."""
    from polyzymd.modifying_proteins import ParameterSet

    if solute.modification_parameters is None:
        return {}
    return {
        site_type: ParameterSet(
            force_field=tuple(chosen.force_field),
            charges=chosen.charges.model_dump(exclude_none=True),
        )
        for site_type, chosen in solute.modification_parameters.as_dict().items()
    }


@dataclass
class ModifiedProteinBuild:
    """A protein after its modifications, with its parameters and records."""

    solute: "ProteinSolute"
    protein: Any
    interchange: "Interchange"
    record: Any
    parameter_sets: dict[str, Any]
    declared: list[dict[str, Any]]
    constraints_from: str
    pdb_path: Path | None = None
    timings: dict[str, float] = field(default_factory=dict)
    #: The config's folder: relative paths in the blocks resolve against it.
    base_dir: Path | None = None

    @classmethod
    def build(
        cls,
        solute: "ProteinSolute",
        constraints_from: str,
        working_dir: Path | None = None,
        base_dir: Path | None = None,
        replicate: int | None = None,
    ) -> "ModifiedProteinBuild":
        """Apply and parameterize the modifications of ``solute``.

        Every change is planned and applied before anything is packed, so a
        site that does not exist or a reaction that does not match stops the
        build here, with the entry it came from.
        """
        from polyzymd import modifying_proteins as mp
        from polyzymd.config.modifications import declared_modifications

        base_dir = Path(base_dir) if base_dir is not None else Path.cwd()
        declared = declared_modifications(
            solute.mutations, solute.modifications, solute.modification_parameters
        )
        timings: dict[str, float] = {}
        t0 = time.perf_counter()
        protein = apply_modifications(solute, base_dir, solute.relax, replicate)
        timings["modify_s"] = time.perf_counter() - t0
        LOGGER.info(
            f"Applied {len(protein.records)} modification(s) to {solute.name} "
            f"({solute.relax} relaxation) in {timings['modify_s']:.1f}s"
        )
        pdb_path = None
        if working_dir is not None:
            pdb_path = Path(working_dir) / MODIFIED_PROTEIN_PDB
            pdb_path.parent.mkdir(parents=True, exist_ok=True)
            protein.save_pdb(pdb_path)
        sets = parameter_sets(solute)
        t0 = time.perf_counter()
        result = mp.parameterize(
            protein,
            sets,
            max_charge_shift=solute.max_charge_shift,
            protein_force_field=tuple(solute.force_field),
            constraints_from=constraints_from,
        )
        # The bare Interchange: GLYCAM06j pairs and masses are written at export.
        interchange, record = result.interchange_without_glycam_pairs, result.record
        timings["parameterize_s"] = time.perf_counter() - t0
        LOGGER.info(
            f"Parameterized {solute.name}: {len(record.sites)} site(s) in "
            f"{timings['parameterize_s']:.1f}s"
        )
        return cls(
            solute=solute,
            protein=protein,
            interchange=interchange,
            record=record,
            parameter_sets=sets,
            declared=declared,
            constraints_from=constraints_from,
            pdb_path=pdb_path,
            timings=timings,
            base_dir=base_dir,
        )

    def molecules(self) -> list["Molecule"]:
        """The protein as OpenFF Molecules, one per connected component in atom order,
        with coordinates and final charges."""
        from openff.toolkit import Molecule
        from openff.units import unit

        from polyzymd.modifying_proteins.parameters import split_molecules

        molecule = Molecule(self.interchange.topology.molecule(0))
        if self.interchange.topology.n_molecules > 1:
            molecule = _joined(self.interchange.topology)
        charges = self.interchange["Electrostatics"].get_charge_array()
        molecule.partial_charges = unit.Quantity(
            [float(q.m_as(unit.elementary_charge)) for q in charges], unit.elementary_charge
        )
        molecule._conformers = [self.interchange.positions.to(unit.angstrom)]
        pieces = split_molecules(molecule)
        for piece in pieces:
            # Residues and chains from each atom's metadata, as Topology.from_pdb
            # gives a protein: without them OpenFF renames atoms to make them unique.
            piece.add_default_hierarchy_schemes()
        return pieces

    @property
    def glycam_sites(self) -> list[Any]:
        """GLYCAM06j sites, whose 1-4 pairs and masses are written at export."""
        return list(self.record.glycam)

    def component(self, name: str, segid: str = "PROT", start: int = 0) -> "ComponentRecord":
        """The protein's component record: one domain for ff14SB, one per site."""
        from polyzymd.modifying_proteins.components import protein_component

        component = protein_component(
            self.protein,
            self.record,
            self.parameter_sets,
            name=name,
            segid=segid,
            start=start,
            base_dir=self.base_dir or ".",
            constraints_from=self.constraints_from,
        )
        # What the config asks for, expanded from the config alone, for the
        # audit to compare with what the build made.
        component.provenance["declared_modifications"] = self.declared
        # The drawn sequence of each grafted chain, for the composition report.
        component.provenance["grafted_chains"] = [
            {
                "segment": record.name,
                "sequence": record.arguments["fragment"]["polymer"]["drawn"],
                "replicate": record.arguments["fragment"]["polymer"].get("replicate"),
            }
            for record in self.protein.records
            if record.kind == "attach"
            and isinstance(record.arguments.get("fragment"), dict)
            and "polymer" in record.arguments["fragment"]
        ]
        return component

    def provenance(self) -> dict[str, Any]:
        """What the build recorded: the changes made, one per site, and how each site was
        parameterized."""
        from polyzymd.builders.polymer_generator import mbuild_provenance
        from polyzymd.modifying_proteins import to_spec

        return {
            "input_pdb": dict(self.protein.source),
            "applied": to_spec(self.protein),
            "declared": self.declared,
            "relax": self.solute.relax,
            "parameterization": self.record.to_dict(),
            "modified_pdb": str(self.pdb_path) if self.pdb_path else None,
            "mbuild": mbuild_provenance(),
            "openff_pablo": _pablo_provenance(),
            "timings_s": {key: round(value, 2) for key, value in self.timings.items()},
        }


def _joined(topology) -> "Molecule":
    """One Molecule holding every molecule of ``topology``, in atom order."""
    from openff.toolkit import Molecule

    joined = Molecule()
    offset = 0
    for molecule in topology.molecules:
        for atom in molecule.atoms:
            joined.add_atom(
                atom.atomic_number,
                atom.formal_charge,
                atom.is_aromatic,
                stereochemistry=atom.stereochemistry,
                name=atom.name,
                metadata=dict(atom.metadata),
            )
        for bond in molecule.bonds:
            joined.add_bond(
                bond.atom1_index + offset,
                bond.atom2_index + offset,
                bond.bond_order,
                bond.is_aromatic,
                stereochemistry=bond.stereochemistry,
            )
        offset += molecule.n_atoms
    return joined


def _pablo_provenance() -> dict[str, Any]:
    import importlib.metadata as metadata
    import json

    try:
        distribution = metadata.distribution("openff-pablo")
    except metadata.PackageNotFoundError:
        return {"package": "openff-pablo", "version": None}
    record: dict[str, Any] = {"package": "openff-pablo", "version": distribution.version}
    direct = distribution.read_text("direct_url.json")
    if direct:
        vcs = json.loads(direct).get("vcs_info", {})
        record.update(commit=vcs.get("commit_id"), requested_revision=vcs.get("requested_revision"))
    return record


def apply_glycam_to_gromacs_export(top_path: Path, sites: list[Any], molecules=None) -> int:
    """Write GLYCAM06j 1-4 pairs and masses into the protein's GROMACS molecule types.

    The protein's molecules come first in ``[ molecules ]``; ``molecules`` are
    their ``[start, stop)`` atom ranges (one molecule spanning every site atom
    when not given). Returns the number of pairs rewritten.
    """
    from polyzymd.modifying_proteins.glycam import apply_to_gromacs_topology

    if not sites:
        return 0
    if molecules is None:
        stop = max(index for site in sites for index in site.atom_indices) + 1
        molecules = [[0, stop]]
    return apply_to_gromacs_topology(top_path, sites, molecules, first_molecule=0)


def apply_modifications(
    solute: "ProteinSolute", base_dir: Path, relax: str, replicate: int | None = None
) -> Any:
    """Load the solute's PDB and apply its blocks; any failure is a ``ValueError``.

    ``replicate`` seeds grafted chains drawn per replicate (see ``config_spec``).

    mBuild raises its own errors (a reaction that matches nothing, a fragment
    it cannot place); they come out as ``ValueError`` with a fix line, as a
    config error does.
    """
    from polyzymd import modifying_proteins as mp

    try:
        protein = mp.load_protein(solute.pdb_path)
        spec = config_spec(solute, replicate)
        if relax == "joint":
            _apply_then_relax_jointly(protein, spec, base_dir)
        else:
            mp.apply_spec(protein, spec, base_dir=base_dir, relax=relax)
    except ValueError:
        raise
    except Exception as exc:  # mBuild's MBuildError, RDKit embedding failures
        raise ValueError(
            f"applying the mutations and modifications of {solute.name} failed: "
            f"{type(exc).__name__}: {exc}\n"
            "fix: check the entry it names (site, atom_name, reaction, fragment); "
            "`polyzymd validate --plans` applies each change to the PDB without building."
        ) from exc
    return protein


def _apply_then_relax_jointly(protein: Any, spec: dict, base_dir: Path) -> None:
    """Place each change with mBuild's torsion search (no minimization), then relax all.

    Rigid placement (what ``apply_spec``'s joint mode does before its one
    relaxation) can drop a dye or a polymer through the protein (sub-angstrom
    heavy-atom pairs) that a minimization then pulls apart by stretching bonds.
    mBuild's torsion search places each fragment and side chain clear of its
    surroundings first; one joint minimization then settles them together,
    as the BG505 build does.
    """
    from polyzymd.modifying_proteins import api
    from polyzymd.modifying_proteins import spec as specs

    mutations, modifications = specs.expand(spec, base_dir)
    for entry in mutations:
        api.mutate(
            protein,
            entry["site"],
            entry["to"],
            stereo=entry.get("stereo"),
            relax=True,
            minimize=False,
        )
    for entry in modifications:
        ((verb, arguments),) = entry.items()
        if verb == "attach":
            arguments = dict(arguments)
            fragment = arguments.pop("fragment")
            arguments.pop("parameters", None)
            api.attach(protein, fragment, relax=True, minimize=False, **arguments)
        else:
            getattr(api, verb)(protein, arguments["site"], arguments["atom_name"])
    api.relax_all(protein)


def plan_rows(solute: "ProteinSolute", base_dir: Path | None = None) -> list[dict[str, Any]]:
    """Return one row per modification the solute's blocks expand to, as the build makes it.

    The blocks are applied to the input structure (without relaxation, so it
    is quick), which checks every site and reaction against it. Each row names
    the modification, its segment, site, what it bonds, the parameter set its
    site gets and the atoms it adds and removes.
    """
    from polyzymd import modifying_proteins as mp
    from polyzymd.utils.mbuild_logs import quiet_mbuild

    # Placement without relaxation leaves close contacts; mBuild warns about
    # them, but a plan is not a structure to simulate, so its log is muted.
    with quiet_mbuild():
        # Plans show replicate 1's grafted chains.
        protein = apply_modifications(solute, base_dir or Path.cwd(), "none", replicate=1)
    sets = solute.modification_parameters.as_dict() if solute.modification_parameters else {}
    owner = {}
    for site in mp.modification_sites(protein):
        for atom in site.atoms:
            owner[atom] = site
    rows = []
    for record in protein.records:
        site = next((owner[a] for a in record.added if a in owner), None)
        if site is None and record.kind != "attach":
            site = next(
                (
                    s
                    for s in set(owner.values())
                    if (record.site.chain_id, record.site.resnum, record.site.icode) in s.residues
                ),
                None,
            )
        chosen = sets.get(site.type) if site is not None else None
        arguments = record.arguments
        rows.append(
            {
                "name": record.name,
                "kind": record.kind,
                "type": record.type,
                "site": str(record.site),
                "residue": f"{record.residue_before}->{record.residue_after}",
                "segment": record.name if record.kind == "attach" else "PROT",
                "fragment": _fragment_label(arguments.get("fragment") or arguments.get("to")),
                "reaction": arguments.get("reaction"),
                "atom_name": arguments.get("atom_name"),
                "parameter_site": site.name if site is not None else None,
                "parameters": (
                    {"force_field": list(chosen.force_field), "charges": chosen.charges.method}
                    if chosen is not None
                    else ("ff14SB" if site is None else None)
                ),
                "atoms_added": len(record.added),
                "atoms_removed": [key.atom_name for key in record.removed],
                "formal_charge_change": record.formal_charge_change,
            }
        )
    return rows


def _fragment_label(source: Any) -> str | None:
    if source is None:
        return None
    if isinstance(source, str):
        return source
    if "smiles" in source:
        return f"{source.get('resname')} {source['smiles']}"
    if "pdb_path" in source:
        return Path(source["pdb_path"]).name
    if "ccd_file" in source:
        return Path(source["ccd_file"]).name
    return source.get("ccd")


def plan_table(rows: list[dict[str, Any]]) -> str:
    """The rows as text, one line per modification."""
    lines = []
    for row in rows:
        what = row["fragment"] or ""
        via = f" via {row['reaction']}" if row["reaction"] else ""
        params = row["parameters"]
        if isinstance(params, dict):
            params = f"{'+'.join(params['force_field'])} ({params['charges']})"
        lines.append(
            f"  {row['segment']:<5} {row['kind']:<11} {row['site']:<10} {row['residue']:<9} "
            f"{what}{via}; site {row['parameter_site'] or '-'}: {params}; "
            f"+{row['atoms_added']} atoms, -{len(row['atoms_removed'])} "
            f"({' '.join(row['atoms_removed']) or 'none'}), "
            f"formal charge {row['formal_charge_change']:+d}"
        )
    return "\n".join(lines)
