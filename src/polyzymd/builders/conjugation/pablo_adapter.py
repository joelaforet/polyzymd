"""Lazy OpenFF Pablo adapter foundation for covalent modification ingestion."""

from __future__ import annotations

import importlib
import importlib.metadata
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.diagnostics import (
    ConjugationDiagnostic,
    DiagnosticCode,
    DiagnosticSeverity,
)
from polyzymd.builders.conjugation.exceptions import (
    ConjugationNotImplementedError,
    PabloIngestionError,
)
from polyzymd.builders.conjugation.metadata import (
    ChainPolicyMetadata,
    ComponentMetadata,
    ComponentRole,
    ConjugationMetadata,
    chain_policy_from_config,
)
from polyzymd.builders.conjugation.structure_inspection import (
    PDBStructureInspection,
    inspect_pdb_structure,
    pdb_atom_records_as_dicts,
)

SUPPORTED_STRUCTURE_SUFFIXES = frozenset({".pdb", ".cif", ".mmcif", ".pdbx"})

STANDARD_PROTEIN_RESIDUES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)
SOLVENT_RESIDUES = frozenset({"HOH", "WAT", "H2O", "SOL", "TIP", "TIP3"})
COMMON_ION_RESIDUES = frozenset(
    {
        "NA",
        "K",
        "CL",
        "CA",
        "MG",
        "ZN",
        "FE",
        "MN",
        "CU",
        "CO",
        "NI",
        "BR",
        "IOD",
    }
)


class PabloAvailability(BaseModel):
    """OpenFF Pablo import availability details."""

    available: bool
    version: str | None = None
    module_path: str | None = None
    error: str | None = None
    warnings: list[str] = Field(default_factory=list)


class PabloStructureCounts(BaseModel):
    """Structure size summary produced during Pablo ingestion."""

    atom_count: int | None = None
    residue_count: int | None = None
    molecule_count: int | None = None
    chain_count: int | None = None
    chain_ids: list[str] = Field(default_factory=list)
    blank_chain_atom_count: int = 0
    blank_chain_residue_count: int = 0
    bond_count: int | None = None
    link_candidate_count: int = 0


class PabloResidueSummary(BaseModel):
    """Residue-level metadata extracted from Pablo or PDB records."""

    residue_id: str
    chain_id: str
    residue_name: str
    residue_number: int | None = None
    residue_index: int | None = None
    atom_count: int = 0
    role: ComponentRole
    is_standard_protein: bool = False
    is_noncanonical: bool = False


class PabloLinkCandidate(BaseModel):
    """Potential covalent junction discovered during ingestion."""

    source: str
    atom_index_1: int | None = None
    atom_name_1: str | None = None
    residue_name_1: str | None = None
    residue_number_1: int | None = None
    chain_id_1: str | None = None
    role_1: ComponentRole | None = None
    atom_index_2: int | None = None
    atom_name_2: str | None = None
    residue_name_2: str | None = None
    residue_number_2: int | None = None
    chain_id_2: str | None = None
    role_2: ComponentRole | None = None
    details: dict[str, Any] = Field(default_factory=dict)


class PabloStructurePreflight(BaseModel):
    """Preflight diagnostics for a structure intended for Pablo ingestion."""

    intended_mode: str
    path: Path
    suffix: str
    pablo: PabloAvailability
    inspection_attempted: bool = False
    inspection_implemented: bool = False
    ingestion_implemented: bool = False
    inspection_summary: PDBStructureInspection | None = None
    warnings: list[str] = Field(default_factory=list)


class PabloIngestionResult(BaseModel):
    """Structured result for one Pablo ingestion attempt."""

    model_config = {"arbitrary_types_allowed": True}

    success: bool
    path: Path
    suffix: str
    pablo: PabloAvailability
    topology: Any | None = Field(default=None, exclude=True)
    counts: PabloStructureCounts = Field(default_factory=PabloStructureCounts)
    residues: list[PabloResidueSummary] = Field(default_factory=list)
    noncanonical_residues: list[PabloResidueSummary] = Field(default_factory=list)
    link_candidates: list[PabloLinkCandidate] = Field(default_factory=list)
    inspection_summary: PDBStructureInspection | None = None
    metadata: ConjugationMetadata = Field(default_factory=ConjugationMetadata)
    diagnostics: list[ConjugationDiagnostic] = Field(default_factory=list)

    def save(self, path: Path | str) -> Path:
        """Save this ingestion result without heavy topology objects.

        Parameters
        ----------
        path : Path or str
            Destination JSON path.

        Returns
        -------
        Path
            Path that was written.
        """
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(self.model_dump(mode="json"), indent=2) + "\n")
        return output_path


class PabloIngestor:
    """Adapter boundary for CCD/Pablo ingestion.

    The class keeps OpenFF Pablo imports inside methods so that importing
    PolyzyMD does not require the full chemistry stack.
    """

    def __init__(self, policy: Any) -> None:
        """Initialize a Pablo ingestion adapter.

        Parameters
        ----------
        policy : Any
            CCD/Pablo policy config from :class:`ConjugationConfig`.
        """
        self._policy = policy

    def ingest_existing(self, pdb_path: Path | str | None) -> Any:
        """Ingest a pre-conjugated structure with OpenFF Pablo.

        Parameters
        ----------
        pdb_path : Path, str, or None
            Path to a pre-conjugated PDB or mmCIF structure.

        Returns
        -------
        Any
            OpenFF topology returned by Pablo.

        Raises
        ------
        ConjugationNotImplementedError
            If Pablo does not return a topology that downstream builders can use.
        PabloIngestionError
            If the input path is missing or invalid.
        """
        result = self.ingest_structure(pdb_path)
        if result.success and result.topology is not None:
            return result.topology

        raise ConjugationNotImplementedError(_not_implemented_message(result))

    def ingest_structure(
        self,
        path: Path | str | None,
        *,
        chain_policy: Any | None = None,
        output_dir: Path | str | None = None,
    ) -> PabloIngestionResult:
        """Attempt real Pablo/CCD ingestion for an existing structure.

        Parameters
        ----------
        path : Path, str, or None
            Structure path to parse. Supported suffixes are PDB and PDBx/mmCIF.
        chain_policy : Any or None, optional
            Chain role policy from the conjugation config, by default ``None``.
        output_dir : Path, str, or None, optional
            Optional directory for a Pablo ingestion sidecar, by default ``None``.

        Returns
        -------
        PabloIngestionResult
            Structured parse result, diagnostics, metadata, and optional topology.

        Raises
        ------
        PabloIngestionError
            If the input path is missing or has an unsupported suffix.
        """
        if path is None:
            raise PabloIngestionError("source_pdb_path is required for Pablo ingestion")

        structure_path = Path(path)
        _validate_structure_path(structure_path)

        chain_metadata = chain_policy_from_config(chain_policy)
        availability = self.probe_available()
        diagnostics = _diagnostics_from_availability(availability)
        inspection = (
            inspect_pdb_structure(structure_path) if _is_pdb_suffix(structure_path) else None
        )
        if inspection is not None:
            diagnostics.extend(_diagnostics_from_inspection(inspection))
        pdb_atom_records = pdb_atom_records_as_dicts(inspection) if inspection is not None else []
        pdb_link_candidates = (
            _pablo_link_candidates_from_inspection(inspection) if inspection is not None else []
        )

        if not availability.available:
            diagnostics.append(
                ConjugationDiagnostic(
                    code=DiagnosticCode.PABLO_INGESTION,
                    severity=DiagnosticSeverity.ERROR,
                    message="OpenFF Pablo is unavailable; structure parsing was not attempted",
                    details={
                        "action": "Run ingestion in the conjugation-py312 pixi environment",
                        "import_error": availability.error,
                    },
                )
            )
            result = _build_result_from_records(
                success=False,
                path=structure_path,
                availability=availability,
                chain_policy=chain_metadata,
                atom_records=pdb_atom_records,
                link_candidates=pdb_link_candidates,
                diagnostics=diagnostics,
                notes=["Pablo was unavailable, so metadata was extracted from PDB records only"],
                inspection=inspection,
            )
            _save_result_sidecar(result, output_dir)
            return result

        try:
            pablo_module = importlib.import_module("openff.pablo")
            topology = pablo_module.topology_from_pdb(
                structure_path,
                residue_library=_build_residue_library(pablo_module, self._policy),
                format=_infer_pablo_format(structure_path),
                use_canonical_names=self._policy_attr("use_canonical_atom_names", False),
            )
        except Exception as exc:  # noqa: BLE001 - third-party errors are normalized to diagnostics
            diagnostics.append(_diagnostic_from_parse_error(exc, structure_path))
            if not pdb_link_candidates:
                diagnostics.append(
                    ConjugationDiagnostic(
                        code=DiagnosticCode.LINK_DISCOVERY,
                        severity=DiagnosticSeverity.WARNING,
                        message=(
                            "Pablo parsing failed before topology bond inspection, and no PDB "
                            "LINK/CONECT cross-role candidates were found"
                        ),
                        details={
                            "action": (
                                "Check whether the input contains LINK/CONECT records for "
                                "nonstandard covalent junctions"
                            )
                        },
                    )
                )
            result = _build_result_from_records(
                success=False,
                path=structure_path,
                availability=availability,
                chain_policy=chain_metadata,
                atom_records=pdb_atom_records,
                link_candidates=pdb_link_candidates,
                diagnostics=diagnostics,
                notes=[
                    "Pablo parsing failed; metadata reflects file-level PDB records when available"
                ],
                inspection=inspection,
            )
            _save_result_sidecar(result, output_dir)
            return result

        topology_atom_records = _atom_records_from_topology(topology)
        topology_link_candidates = _topology_link_candidates(topology, chain_metadata)
        link_candidates = [*topology_link_candidates, *pdb_link_candidates]
        if not topology_link_candidates:
            diagnostics.append(
                ConjugationDiagnostic(
                    code=DiagnosticCode.LINK_DISCOVERY,
                    severity=DiagnosticSeverity.WARNING,
                    message=(
                        "Pablo parsed the structure, but no cross-role topology bonds were "
                        "identified as covalent link candidates"
                    ),
                    details={
                        "pdb_link_or_conect_candidates": len(pdb_link_candidates),
                        "action": (
                            "Review component chain policy and PDB LINK/CONECT records if a "
                            "protein-moiety junction is expected"
                        ),
                    },
                )
            )
        diagnostics.append(
            ConjugationDiagnostic(
                code=DiagnosticCode.PABLO_INGESTION,
                severity=DiagnosticSeverity.INFO,
                message="OpenFF Pablo parsed the structure into an OpenFF topology",
                details={"topology_type": type(topology).__name__},
            )
        )
        result = _build_result_from_records(
            success=True,
            path=structure_path,
            availability=availability,
            chain_policy=chain_metadata,
            atom_records=topology_atom_records,
            link_candidates=link_candidates,
            diagnostics=diagnostics,
            topology=topology,
            molecule_count=_safe_int_attr(topology, "n_molecules"),
            bond_count=_safe_int_attr(topology, "n_bonds"),
            notes=["Pablo topology is available for downstream parameterization phases"],
            inspection=inspection,
        )
        _save_result_sidecar(result, output_dir)
        return result

    def check_available(self) -> PabloAvailability:
        """Import OpenFF Pablo lazily and return version/path details.

        Returns
        -------
        PabloAvailability
            Successful Pablo import details.

        Raises
        ------
        PabloIngestionError
            If OpenFF Pablo is not importable in the active environment.
        """
        availability = self.probe_available()
        if not availability.available:
            raise PabloIngestionError(
                "OpenFF Pablo is not importable in this environment. Use the "
                "'conjugation-py312' pixi environment for conjugation ingestion preflight. "
                f"Original import error: {availability.error}"
            )
        return availability

    def probe_available(self) -> PabloAvailability:
        """Probe OpenFF Pablo availability without raising on import failure.

        Returns
        -------
        PabloAvailability
            Pablo import details or a structured import error.
        """
        try:
            pablo_module = importlib.import_module("openff.pablo")
        except ImportError as exc:
            return PabloAvailability(
                available=False,
                error=str(exc),
                warnings=[
                    "OpenFF Pablo is unavailable; install/use the conjugation-py312 pixi "
                    "environment before production ingestion"
                ],
            )

        version = _detect_pablo_version(pablo_module)
        module_path = getattr(pablo_module, "__file__", None)
        return PabloAvailability(
            available=True,
            version=version,
            module_path=str(module_path) if module_path is not None else None,
            warnings=["Pablo import succeeded; structure ingestion can be attempted"],
        )

    def preflight_structure(self, path: Path | str) -> PabloStructurePreflight:
        """Validate a structure path and report Pablo readiness.

        Parameters
        ----------
        path : Path or str
            Structure path intended for Pablo ingestion.

        Returns
        -------
        PabloStructurePreflight
            File validation and Pablo availability diagnostics.

        Raises
        ------
        PabloIngestionError
            If the path is missing or has an unsupported suffix.
        """
        structure_path = Path(path)
        _validate_structure_path(structure_path)

        availability = self.probe_available()
        inspection = (
            inspect_pdb_structure(structure_path) if _is_pdb_suffix(structure_path) else None
        )
        warnings = [
            "Pablo structure ingestion is available through ingest_structure()",
            *availability.warnings,
        ]
        if inspection is not None:
            warnings.extend(inspection.convention_warnings)
            warnings.extend(inspection.compatibility_warnings)
        return PabloStructurePreflight(
            intended_mode="ingest_existing",
            path=structure_path,
            suffix=structure_path.suffix.lower(),
            pablo=availability,
            inspection_attempted=inspection is not None,
            inspection_implemented=True,
            ingestion_implemented=True,
            inspection_summary=inspection,
            warnings=warnings,
        )

    def _policy_attr(self, name: str, default: Any) -> Any:
        """Return an attribute from the configured Pablo policy."""
        return _policy_attr_from(self._policy, name, default)


def _validate_structure_path(path: Path) -> None:
    """Validate a structure path for Pablo ingestion.

    Parameters
    ----------
    path : Path
        Structure path to validate.

    Raises
    ------
    PabloIngestionError
        If the path does not exist, is not a file, or has an unsupported suffix.
    """
    if not path.exists():
        raise PabloIngestionError(f"Structure file does not exist: {path}")
    if not path.is_file():
        raise PabloIngestionError(f"Structure path is not a file: {path}")
    if path.suffix.lower() not in SUPPORTED_STRUCTURE_SUFFIXES:
        supported = ", ".join(sorted(SUPPORTED_STRUCTURE_SUFFIXES))
        raise PabloIngestionError(
            f"Unsupported structure suffix '{path.suffix}'. Expected one of: {supported}"
        )


def _detect_pablo_version(pablo_module: Any) -> str | None:
    """Return the OpenFF Pablo version when available.

    Parameters
    ----------
    pablo_module : Any
        Imported Pablo module.

    Returns
    -------
    str or None
        Package or module version.
    """
    try:
        return importlib.metadata.version("openff-pablo")
    except importlib.metadata.PackageNotFoundError:
        return getattr(pablo_module, "__version__", None)


def _build_residue_library(pablo_module: Any, policy: Any | None) -> Any:
    """Build the residue library used by Pablo parsing.

    Parameters
    ----------
    pablo_module : Any
        Imported Pablo module.
    policy : Any or None
        CCD/Pablo policy config.

    Returns
    -------
    Any
        Pablo residue library mapping.
    """
    residue_library = pablo_module.STD_CCD_CACHE
    lookup_policy = _policy_attr_from(policy, "lookup_policy", None)
    lookup_value = getattr(lookup_policy, "value", lookup_policy)
    if lookup_value == "auto_download" and hasattr(residue_library, "with_"):
        residue_library = residue_library.with_({})
        residue_library.auto_download = True
    return residue_library


def _policy_attr_from(policy: Any | None, name: str, default: Any) -> Any:
    """Return a policy attribute without depending on config model imports."""
    if policy is None:
        return default
    return getattr(policy, name, default)


def _infer_pablo_format(path: Path) -> str:
    """Infer Pablo's explicit format string from a path suffix."""
    if path.suffix.lower() in {".cif", ".mmcif", ".pdbx"}:
        return "CIF"
    return "PDB"


def _is_pdb_suffix(path: Path) -> bool:
    """Return whether a path should be parsed as fixed-width PDB records."""
    return path.suffix.lower() == ".pdb"


def _diagnostics_from_availability(availability: PabloAvailability) -> list[ConjugationDiagnostic]:
    """Create diagnostics from a Pablo availability result."""
    severity = DiagnosticSeverity.INFO if availability.available else DiagnosticSeverity.WARNING
    return [
        ConjugationDiagnostic(
            code=DiagnosticCode.PABLO_ADAPTER,
            severity=severity,
            message="OpenFF Pablo availability checked for structure ingestion",
            details=availability.model_dump(mode="json"),
        )
    ]


def _diagnostics_from_inspection(
    inspection: PDBStructureInspection,
) -> list[ConjugationDiagnostic]:
    """Create diagnostics from pure-Python PDB inspection."""
    diagnostics = [
        ConjugationDiagnostic(
            code=DiagnosticCode.PDB_STRUCTURE_INSPECTION,
            severity=DiagnosticSeverity.INFO,
            message="PDB structure inspection completed for Pablo preflight",
            details={
                "atom_count": inspection.atom_count,
                "residue_count": inspection.residue_count,
                "chain_ids": inspection.chain_ids,
                "blank_chain_atom_count": inspection.blank_chain_atom_count,
                "blank_chain_residue_count": inspection.blank_chain_residue_count,
                "noncanonical_residue_count": len(inspection.noncanonical_residue_candidates),
                "polymer_ptm_candidate_count": len(inspection.polymer_ptm_candidates),
                "covalent_attachment_candidate_count": len(
                    inspection.covalent_attachment_candidates
                ),
                "ssbond_count": inspection.ssbond_count,
                "residue_name_counts": inspection.residue_name_counts,
            },
        )
    ]
    warnings = [*inspection.convention_warnings, *inspection.compatibility_warnings]
    if warnings:
        diagnostics.append(
            ConjugationDiagnostic(
                code=DiagnosticCode.PDB_STRUCTURE_INSPECTION,
                severity=DiagnosticSeverity.WARNING,
                message="PDB/Pablo compatibility warnings were found",
                details={"warnings": warnings},
            )
        )
    return diagnostics


def _pablo_link_candidates_from_inspection(
    inspection: PDBStructureInspection,
) -> list[PabloLinkCandidate]:
    """Convert PDB inspection attachment evidence to Pablo adapter candidates."""
    candidates: list[PabloLinkCandidate] = []
    for candidate in inspection.covalent_attachment_candidates:
        role_2 = (
            ComponentRole.POLYMER
            if candidate.candidate_category == "polymer_ptm"
            else ComponentRole.MOIETY
        )
        candidates.append(
            PabloLinkCandidate(
                source=candidate.source,
                atom_name_1=candidate.protein_atom_name,
                residue_name_1=candidate.protein_residue_name,
                residue_number_1=_residue_number_from_id(candidate.protein_residue_id),
                chain_id_1=candidate.protein_chain_id,
                role_1=ComponentRole.PROTEIN,
                atom_name_2=candidate.candidate_atom_name,
                residue_name_2=candidate.candidate_residue_name,
                residue_number_2=_residue_number_from_id(candidate.candidate_residue_id),
                chain_id_2=candidate.candidate_chain_id,
                role_2=role_2,
                details={
                    **candidate.details,
                    "line_number": candidate.line_number,
                    "distance_angstrom": candidate.distance_angstrom,
                    "candidate_category": candidate.candidate_category,
                },
            )
        )
    return candidates


def _diagnostic_from_parse_error(exc: Exception, path: Path) -> ConjugationDiagnostic:
    """Convert a Pablo parse exception into an actionable diagnostic."""
    return ConjugationDiagnostic(
        code=DiagnosticCode.PABLO_INGESTION,
        severity=DiagnosticSeverity.ERROR,
        message="OpenFF Pablo could not parse the structure into a topology",
        details={
            "path": str(path),
            "exception_type": type(exc).__name__,
            "error": str(exc),
            "action": (
                "Ensure all atoms and hydrogens are present, residue and atom names match CCD "
                "definitions, and LINK/CONECT records describe nonstandard covalent bonds. "
                "Custom residue definition wiring is deferred to a later conjugation phase."
            ),
        },
    )


def _read_pdb_atom_records(path: Path) -> list[dict[str, Any]]:
    """Read lightweight ATOM/HETATM metadata from a PDB file."""
    records: list[dict[str, Any]] = []
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), start=1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue_number = _parse_int(line[22:26].strip())
        records.append(
            {
                "atom_index": len(records),
                "atom_serial": line[6:11].strip() or None,
                "atom_name": line[12:16].strip() or None,
                "residue_name": line[17:20].strip() or "",
                "chain_id": line[21:22].strip() or "",
                "residue_number": residue_number,
                "res_seq": line[22:26].strip() or None,
                "residue_index": None,
                "insertion_code": line[26:27].strip() or None,
                "record_name": line[:6].strip(),
                "line_number": line_number,
            }
        )
    return records


def _read_pdb_link_candidates(
    path: Path,
    chain_policy: ChainPolicyMetadata,
    atom_records: list[dict[str, Any]],
) -> list[PabloLinkCandidate]:
    """Extract LINK and cross-role CONECT records from a PDB file."""
    if not _is_pdb_suffix(path):
        return []

    serial_records = {
        str(record.get("atom_serial")): record
        for record in atom_records
        if record.get("atom_serial") is not None
    }
    candidates: list[PabloLinkCandidate] = []
    seen_conect_pairs: set[tuple[str, str]] = set()
    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), start=1):
        if line.startswith("LINK"):
            candidates.append(_link_candidate_from_link_record(line, line_number, chain_policy))
        elif line.startswith("CONECT"):
            source_serial = line[6:11].strip()
            for target_serial in _parse_conect_targets(line):
                pair = tuple(sorted((source_serial, target_serial)))
                if pair in seen_conect_pairs:
                    continue
                seen_conect_pairs.add(pair)
                source_record = serial_records.get(source_serial)
                target_record = serial_records.get(target_serial)
                if source_record is None or target_record is None:
                    continue
                candidate = _link_candidate_from_atom_records(
                    source="CONECT",
                    record_1=source_record,
                    record_2=target_record,
                    chain_policy=chain_policy,
                    details={"line_number": line_number, "atom_serials": list(pair)},
                )
                if _is_cross_role_candidate(candidate):
                    candidates.append(candidate)
    return candidates


def _parse_conect_targets(line: str) -> list[str]:
    """Parse bonded atom serials from a CONECT record."""
    targets: list[str] = []
    for start in range(11, len(line), 5):
        target = line[start : start + 5].strip()
        if target:
            targets.append(target)
    return targets


def _link_candidate_from_link_record(
    line: str,
    line_number: int,
    chain_policy: ChainPolicyMetadata,
) -> PabloLinkCandidate:
    """Build a link candidate from a fixed-width LINK record."""
    residue_name_1 = line[17:20].strip() or None
    chain_id_1 = line[21:22].strip() or None
    residue_name_2 = line[47:50].strip() or None
    chain_id_2 = line[51:52].strip() or None
    role_1 = _role_for_record(chain_id_1 or "", residue_name_1 or "", chain_policy)
    role_2 = _role_for_record(chain_id_2 or "", residue_name_2 or "", chain_policy)
    return PabloLinkCandidate(
        source="LINK",
        atom_name_1=line[12:16].strip() or None,
        residue_name_1=residue_name_1,
        residue_number_1=_parse_int(line[22:26].strip()),
        chain_id_1=chain_id_1,
        role_1=role_1,
        atom_name_2=line[42:46].strip() or None,
        residue_name_2=residue_name_2,
        residue_number_2=_parse_int(line[52:56].strip()),
        chain_id_2=chain_id_2,
        role_2=role_2,
        details={"line_number": line_number, "raw_record": line.rstrip()},
    )


def _atom_records_from_topology(topology: Any) -> list[dict[str, Any]]:
    """Extract atom and residue metadata from an OpenFF topology."""
    records: list[dict[str, Any]] = []
    for atom_index, atom in enumerate(topology.atoms):
        metadata = getattr(atom, "metadata", {}) or {}
        records.append(
            {
                "atom_index": atom_index,
                "atom_serial": metadata.get("atom_serial"),
                "atom_name": getattr(atom, "name", None) or metadata.get("used_synonym"),
                "residue_name": metadata.get("residue_name", "") or "",
                "chain_id": metadata.get("chain_id", "") or "",
                "residue_number": _parse_int(metadata.get("residue_number")),
                "res_seq": metadata.get("res_seq"),
                "residue_index": _parse_int(metadata.get("residue_index")),
                "insertion_code": metadata.get("insertion_code"),
                "pdb_index": metadata.get("pdb_index"),
            }
        )
    return records


def _topology_link_candidates(
    topology: Any,
    chain_policy: ChainPolicyMetadata,
) -> list[PabloLinkCandidate]:
    """Extract cross-role bond candidates from a Pablo/OpenFF topology."""
    atoms = list(topology.atoms)
    atom_records = _atom_records_from_topology(topology)
    atom_indices = {id(atom): index for index, atom in enumerate(atoms)}
    candidates: list[PabloLinkCandidate] = []
    for bond in topology.bonds:
        atom_index_1 = atom_indices.get(id(bond.atom1))
        atom_index_2 = atom_indices.get(id(bond.atom2))
        if atom_index_1 is None or atom_index_2 is None:
            continue
        candidate = _link_candidate_from_atom_records(
            source="PabloTopologyBond",
            record_1=atom_records[atom_index_1],
            record_2=atom_records[atom_index_2],
            chain_policy=chain_policy,
            details={"bond_type": type(bond).__name__},
        )
        if _is_cross_role_candidate(candidate):
            candidates.append(candidate)
    return candidates


def _link_candidate_from_atom_records(
    *,
    source: str,
    record_1: dict[str, Any],
    record_2: dict[str, Any],
    chain_policy: ChainPolicyMetadata,
    details: dict[str, Any] | None = None,
) -> PabloLinkCandidate:
    """Build a link candidate from two atom metadata records."""
    role_1 = _role_for_record(
        str(record_1.get("chain_id") or ""),
        str(record_1.get("residue_name") or ""),
        chain_policy,
    )
    role_2 = _role_for_record(
        str(record_2.get("chain_id") or ""),
        str(record_2.get("residue_name") or ""),
        chain_policy,
    )
    return PabloLinkCandidate(
        source=source,
        atom_index_1=_parse_int(record_1.get("atom_index")),
        atom_name_1=record_1.get("atom_name"),
        residue_name_1=record_1.get("residue_name"),
        residue_number_1=_parse_int(record_1.get("residue_number")),
        chain_id_1=record_1.get("chain_id"),
        role_1=role_1,
        atom_index_2=_parse_int(record_2.get("atom_index")),
        atom_name_2=record_2.get("atom_name"),
        residue_name_2=record_2.get("residue_name"),
        residue_number_2=_parse_int(record_2.get("residue_number")),
        chain_id_2=record_2.get("chain_id"),
        role_2=role_2,
        details=details or {},
    )


def _is_cross_role_candidate(candidate: PabloLinkCandidate) -> bool:
    """Return whether a link candidate crosses component roles."""
    if candidate.role_1 is None or candidate.role_2 is None:
        return False
    if candidate.role_1 == candidate.role_2:
        return False
    covalent_roles = {ComponentRole.PROTEIN, ComponentRole.MOIETY, ComponentRole.POLYMER}
    return candidate.role_1 in covalent_roles or candidate.role_2 in covalent_roles


def _build_result_from_records(
    *,
    success: bool,
    path: Path,
    availability: PabloAvailability,
    chain_policy: ChainPolicyMetadata,
    atom_records: list[dict[str, Any]],
    link_candidates: list[PabloLinkCandidate],
    diagnostics: list[ConjugationDiagnostic],
    notes: list[str],
    topology: Any | None = None,
    molecule_count: int | None = None,
    bond_count: int | None = None,
    inspection: PDBStructureInspection | None = None,
) -> PabloIngestionResult:
    """Build a complete ingestion result from atom-level records."""
    residues = _summarize_residues(atom_records, chain_policy)
    noncanonical = [residue for residue in residues if residue.is_noncanonical]
    components = _build_components(atom_records, chain_policy)
    chain_ids = sorted({str(record.get("chain_id") or "") for record in atom_records})
    counts = PabloStructureCounts(
        atom_count=len(atom_records) if atom_records else None,
        residue_count=len(residues) if residues else None,
        molecule_count=molecule_count,
        chain_count=len(chain_ids) if atom_records else None,
        chain_ids=chain_ids,
        blank_chain_atom_count=sum(1 for record in atom_records if not record.get("chain_id")),
        blank_chain_residue_count=inspection.blank_chain_residue_count
        if inspection is not None
        else 0,
        bond_count=bond_count,
        link_candidate_count=len(link_candidates),
    )
    metadata = ConjugationMetadata(
        enabled=True,
        mode="ingest_existing",
        chain_policy=chain_policy,
        components=components,
        notes=notes,
    )
    diagnostics.append(
        ConjugationDiagnostic(
            code=DiagnosticCode.COMPONENT_METADATA,
            severity=DiagnosticSeverity.INFO,
            message="Component metadata extracted for Pablo ingestion",
            details={
                "atom_count": counts.atom_count,
                "residue_count": counts.residue_count,
                "chain_count": counts.chain_count,
                "noncanonical_residue_count": len(noncanonical),
            },
        )
    )
    return PabloIngestionResult(
        success=success,
        path=path,
        suffix=path.suffix.lower(),
        pablo=availability,
        topology=topology,
        counts=counts,
        residues=residues,
        noncanonical_residues=noncanonical,
        link_candidates=link_candidates,
        inspection_summary=inspection,
        metadata=metadata,
        diagnostics=diagnostics,
    )


def _summarize_residues(
    atom_records: list[dict[str, Any]],
    chain_policy: ChainPolicyMetadata,
) -> list[PabloResidueSummary]:
    """Summarize atom-level records into residue metadata."""
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for record in atom_records:
        key = (
            record.get("chain_id", ""),
            record.get("residue_name", ""),
            record.get("residue_number"),
            record.get("res_seq"),
            record.get("residue_index"),
            record.get("insertion_code"),
        )
        grouped[key].append(record)

    residues: list[PabloResidueSummary] = []
    for fallback_index, (key, records) in enumerate(grouped.items()):
        chain_id, residue_name, residue_number, res_seq, residue_index, insertion_code = key
        role = _role_for_record(chain_id, residue_name, chain_policy)
        residue_id = _format_residue_id(
            chain_id=chain_id,
            residue_name=residue_name,
            residue_number=residue_number,
            res_seq=res_seq,
            residue_index=residue_index if residue_index is not None else fallback_index,
            insertion_code=insertion_code,
        )
        residues.append(
            PabloResidueSummary(
                residue_id=residue_id,
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=_parse_int(residue_number),
                residue_index=_parse_int(residue_index),
                atom_count=len(records),
                role=role,
                is_standard_protein=residue_name in STANDARD_PROTEIN_RESIDUES,
                is_noncanonical=_is_noncanonical_residue(residue_name),
            )
        )
    return residues


def _build_components(
    atom_records: list[dict[str, Any]],
    chain_policy: ChainPolicyMetadata,
) -> list[ComponentMetadata]:
    """Build component metadata from atom-level records."""
    grouped: dict[tuple[ComponentRole, str], list[dict[str, Any]]] = defaultdict(list)
    for record in atom_records:
        chain_id = str(record.get("chain_id") or "")
        residue_name = str(record.get("residue_name") or "")
        role = _role_for_record(chain_id, residue_name, chain_policy)
        grouped[(role, chain_id)].append(record)

    components: list[ComponentMetadata] = []
    for (role, chain_id), records in sorted(
        grouped.items(),
        key=lambda item: (item[0][1], item[0][0]),
    ):
        residue_numbers = sorted(
            {
                number
                for number in (_parse_int(record.get("residue_number")) for record in records)
                if number is not None
            }
        )
        residue_names = sorted({str(record.get("residue_name") or "") for record in records})
        noncanonical_names = [name for name in residue_names if _is_noncanonical_residue(name)]
        component_chain = chain_id or "blank"
        components.append(
            ComponentMetadata(
                component_id=f"{role.value}:{component_chain}",
                role=role,
                chain_id=chain_id,
                atom_indices=[_parse_int(record.get("atom_index")) or 0 for record in records],
                residue_numbers=residue_numbers,
                source="Pablo" if any("pdb_index" in record for record in records) else "PDB",
                details={
                    "atom_count": len(records),
                    "residue_count": len(
                        {
                            (
                                record.get("residue_name"),
                                record.get("residue_number"),
                                record.get("res_seq"),
                                record.get("residue_index"),
                            )
                            for record in records
                        }
                    ),
                    "residue_names": residue_names,
                    "noncanonical_residue_names": noncanonical_names,
                },
            )
        )
    return components


def _role_for_record(
    chain_id: str,
    residue_name: str,
    chain_policy: ChainPolicyMetadata,
) -> ComponentRole:
    """Assign a component role from chain policy and residue identity."""
    normalized_chain = (chain_id or "").upper()
    normalized_residue = (residue_name or "").upper()
    if normalized_chain == chain_policy.protein_chain:
        return ComponentRole.PROTEIN
    if normalized_chain == chain_policy.substrate_chain:
        return ComponentRole.SUBSTRATE
    if normalized_chain == chain_policy.moiety_chain:
        return ComponentRole.MOIETY
    if normalized_residue in SOLVENT_RESIDUES or normalized_residue in COMMON_ION_RESIDUES:
        return ComponentRole.SOLVENT
    if _is_noncanonical_residue(normalized_residue):
        return ComponentRole.MOIETY
    if normalized_chain >= chain_policy.solvent_start_chain:
        return ComponentRole.SOLVENT
    return ComponentRole.PROTEIN


def _is_noncanonical_residue(residue_name: str) -> bool:
    """Return whether a residue name is noncanonical for protein ingestion."""
    normalized = (residue_name or "").upper()
    if not normalized:
        return False
    if normalized in STANDARD_PROTEIN_RESIDUES:
        return False
    if normalized in SOLVENT_RESIDUES or normalized in COMMON_ION_RESIDUES:
        return False
    return True


def _format_residue_id(
    *,
    chain_id: str,
    residue_name: str,
    residue_number: Any,
    res_seq: Any,
    residue_index: Any,
    insertion_code: Any,
) -> str:
    """Format a stable residue identifier for diagnostics."""
    number = residue_number if residue_number is not None else res_seq
    if number in (None, ""):
        number = f"idx{residue_index}"
    insertion = insertion_code or ""
    return f"{chain_id or '_'}:{residue_name}{number}{insertion}"


def _residue_number_from_id(residue_id: str) -> int | None:
    """Extract a residue number from a formatted diagnostics residue ID."""
    residue_label = residue_id.split(":", maxsplit=1)[-1]
    digits = "".join(character for character in residue_label if character.isdigit())
    return _parse_int(digits)


def _parse_int(value: Any) -> int | None:
    """Parse an integer-like value and return ``None`` on failure."""
    if value is None:
        return None
    if isinstance(value, int):
        return value
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _safe_int_attr(obj: Any, attr: str) -> int | None:
    """Return an integer attribute when present."""
    return _parse_int(getattr(obj, attr, None))


def _save_result_sidecar(result: PabloIngestionResult, output_dir: Path | str | None) -> None:
    """Save the optional Pablo ingestion sidecar."""
    if output_dir is None:
        return
    result.save(Path(output_dir) / "pablo_ingestion_result.json")


def _not_implemented_message(result: PabloIngestionResult) -> str:
    """Build a clear downstream-not-implemented message from an ingestion result."""
    error_details = [
        diagnostic.details.get("error")
        for diagnostic in result.diagnostics
        if diagnostic.severity == DiagnosticSeverity.ERROR and diagnostic.details.get("error")
    ]
    suffix = f" Last Pablo error: {error_details[0]}" if error_details else ""
    return (
        "OpenFF Pablo ingestion did not produce a usable topology for downstream "
        "conjugation parameterization. Review conjugation diagnostics and the "
        "pablo_ingestion_result.json sidecar for actionable parser details."
        f"{suffix}"
    )
