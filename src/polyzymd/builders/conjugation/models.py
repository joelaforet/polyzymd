"""Public models for conjugate construction."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field, model_validator

from polyzymd.config.schema import (
    ConjugationAttachmentConfig,
    ConjugationCcdPabloPolicyConfig,
    ConjugationChainPolicyConfig,
)


class ConjugateBuildRequest(BaseModel):
    """Lightweight public request shell for conjugate construction."""

    model_config = {"arbitrary_types_allowed": True}

    config: Any | None = Field(
        default=None,
        description="In-memory SimulationConfig or compatible config object.",
    )
    config_path: Path | None = Field(
        default=None,
        description="Path to a YAML config consumed by the conjugation workflow.",
    )
    output_dir: Path | None = Field(
        default=None,
        description="Output directory for generated conjugation artifacts.",
    )
    free_polymer_seed: int | None = Field(
        default=None,
        description="Optional seed forwarded to free-polymer placement.",
    )
    protein_pdb_path: Path | None = Field(
        default=None,
        description="Cleaned unmodified protein PDB for direct public conjugation requests.",
    )
    attachments: tuple[ConjugationAttachmentConfig, ...] = Field(
        default_factory=tuple,
        description="Covalent modifications to apply to the cleaned protein.",
    )
    ccd_pablo: ConjugationCcdPabloPolicyConfig | None = Field(
        default=None,
        description="Optional Pablo/CCD policy for direct conjugation requests.",
    )
    chain_policy: ConjugationChainPolicyConfig | None = Field(
        default=None,
        description="Optional chain assignment policy for direct conjugation requests.",
    )

    @model_validator(mode="after")
    def validate_supported_sources(self) -> "ConjugateBuildRequest":
        """Require either a config source or the direct protein/modification source."""
        has_config = self.config is not None or self.config_path is not None
        has_direct = self.protein_pdb_path is not None or bool(self.attachments)
        if has_config and has_direct:
            raise ValueError(
                "ConjugateBuildRequest accepts either config/config_path or "
                "protein_pdb_path with attachments, not both"
            )
        if not has_config and has_direct:
            if self.protein_pdb_path is None:
                raise ValueError("protein_pdb_path is required when attachments are supplied")
            if not self.attachments:
                raise ValueError("attachments are required when protein_pdb_path is supplied")
        return self


class ConjugationResult(BaseModel):
    """Public conjugation orchestration result.

    The result carries lightweight artifact paths for serialization and keeps
    heavyweight build objects in excluded in-memory fields for follow-on export.
    """

    model_config = {"arbitrary_types_allowed": True}

    status: str = "completed"
    output_dir: Path | None = None
    config_path: Path | None = None
    crosslinked_conjugate_pdb_path: Path | None = None
    relaxed_conjugate_pdb_path: Path | None = None
    solvated_pdb_path: Path | None = None
    workflow_json_path: Path | None = None
    final_interchange_created: bool | None = None
    artifact_paths: dict[str, Path] = Field(default_factory=dict)
    generated_sequence: str | None = None
    reactive_sequence_index: int | None = None
    reactive_residue_selector: dict[str, int | str] | None = None
    generated_sequences: tuple[str, ...] = Field(default_factory=tuple)
    reactive_sequence_indices: tuple[int, ...] = Field(default_factory=tuple)
    reactive_residue_selectors: tuple[dict[str, int | str], ...] = Field(default_factory=tuple)
    conjugate_generation: Any | None = Field(default=None, exclude=True)
    conjugate_generations: tuple[Any, ...] = Field(default_factory=tuple, exclude=True)
    construction: Any | None = Field(default=None, exclude=True)
    attachment_specs: tuple[Any, ...] = Field(default_factory=tuple, exclude=True)
    protein_canonicalization: Any | None = Field(default=None, exclude=True)
    modifier: Any | None = Field(default=None, exclude=True)
    modifiers: tuple[Any, ...] = Field(default_factory=tuple, exclude=True)
    final_interchange: Any | None = Field(default=None, exclude=True)
    exact_export_bundle: Any | None = Field(default=None, exclude=True)
    system_builder: Any | None = Field(default=None, exclude=True)
    relaxed_conjugate_topology: Any = Field(default=None, exclude=True)
    solvated_topology: Any = Field(default=None, exclude=True)

    @model_validator(mode="after")
    def populate_artifact_paths(self) -> "ConjugationResult":
        """Populate lightweight path fields from retained workflow objects."""
        construction = self.construction
        if self.crosslinked_conjugate_pdb_path is None:
            self.crosslinked_conjugate_pdb_path = _optional_path(
                getattr(construction, "crosslinked_pdb_path", None)
            )
        paths = {
            "crosslinked_conjugate_pdb": self.crosslinked_conjugate_pdb_path,
            "relaxed_conjugate_pdb": self.relaxed_conjugate_pdb_path,
            "solvated_pdb": self.solvated_pdb_path,
            "workflow_json": self.workflow_json_path,
            "conjugate_validation_report": _optional_path(
                getattr(construction, "validation_report_path", None)
            ),
        }
        self.artifact_paths = {
            **{name: path for name, path in paths.items() if path is not None},
            **self.artifact_paths,
        }
        _restore_legacy_pdb_fragment_alias(self.artifact_paths)
        return self

    @classmethod
    def from_workflow_result(
        cls,
        workflow_result: Any,
        *,
        config_path: Path | str | None = None,
    ) -> "ConjugationResult":
        """Create a public result from a workflow-compatible result object.

        Parameters
        ----------
        workflow_result : Any
            Result object produced by a conjugation workflow.
        config_path : Path or str, optional
            Config path used for path-driven builds.

        Returns
        -------
        ConjugationResult
            Public build result with artifact paths and excluded heavy objects.
        """
        if isinstance(workflow_result, cls):
            if config_path is None:
                return workflow_result
            return workflow_result.model_copy(update={"config_path": _optional_path(config_path)})

        construction = getattr(workflow_result, "construction", None)
        crosslinked_path = _optional_path(getattr(construction, "crosslinked_pdb_path", None))
        relaxed_path = _optional_path(getattr(workflow_result, "relaxed_conjugate_pdb_path", None))
        solvated_path = _optional_path(getattr(workflow_result, "solvated_pdb_path", None))
        workflow_path = _optional_path(getattr(workflow_result, "workflow_json_path", None))
        validation_path = _optional_path(getattr(construction, "validation_report_path", None))
        artifact_paths = {
            name: path
            for name, path in {
                "crosslinked_conjugate_pdb": crosslinked_path,
                "relaxed_conjugate_pdb": relaxed_path,
                "solvated_pdb": solvated_path,
                "workflow_json": workflow_path,
                "conjugate_validation_report": validation_path,
            }.items()
            if path is not None
        }
        return cls(
            status="completed",
            output_dir=_optional_path(getattr(workflow_result, "output_dir", None)),
            config_path=_optional_path(config_path),
            crosslinked_conjugate_pdb_path=crosslinked_path,
            relaxed_conjugate_pdb_path=relaxed_path,
            solvated_pdb_path=solvated_path,
            workflow_json_path=workflow_path,
            final_interchange_created=getattr(workflow_result, "final_interchange_created", None),
            artifact_paths=artifact_paths,
            generated_sequence=getattr(workflow_result, "generated_sequence", None),
            reactive_sequence_index=getattr(workflow_result, "reactive_sequence_index", None),
            reactive_residue_selector=getattr(workflow_result, "reactive_residue_selector", None),
            generated_sequences=tuple(getattr(workflow_result, "generated_sequences", ()) or ()),
            reactive_sequence_indices=tuple(
                getattr(workflow_result, "reactive_sequence_indices", ()) or ()
            ),
            reactive_residue_selectors=tuple(
                getattr(workflow_result, "reactive_residue_selectors", ()) or ()
            ),
            conjugate_generation=getattr(workflow_result, "conjugate_generation", None),
            conjugate_generations=tuple(
                getattr(workflow_result, "conjugate_generations", ()) or ()
            ),
            construction=construction,
            attachment_specs=tuple(getattr(workflow_result, "attachment_specs", ()) or ()),
            protein_canonicalization=getattr(workflow_result, "protein_canonicalization", None),
            modifier=getattr(workflow_result, "modifier", None),
            modifiers=tuple(getattr(workflow_result, "modifiers", ()) or ()),
            relaxed_conjugate_topology=getattr(workflow_result, "relaxed_conjugate_topology", None),
            solvated_topology=getattr(workflow_result, "solvated_topology", None),
            final_interchange=getattr(workflow_result, "final_interchange", None),
            exact_export_bundle=getattr(workflow_result, "exact_export_bundle", None),
            system_builder=getattr(workflow_result, "system_builder", None),
        )

    def save(self, path: Path | str) -> Path:
        """Write a JSON sidecar excluding heavy in-memory build objects.

        Parameters
        ----------
        path : Path or str
            Destination JSON path.

        Returns
        -------
        Path
            Written destination path.
        """
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(self.model_dump(mode="json"), indent=2) + "\n")
        return output_path

    def require_final_interchange(self) -> Any:
        """Return the final Interchange object or fail with an actionable error.

        Returns
        -------
        Any
            In-memory final Interchange object.

        Raises
        ------
        RuntimeError
            If the workflow did not create or retain a final Interchange.
        """
        if self.final_interchange is None:
            if self.exact_export_bundle is not None:
                return self.exact_export_bundle
            raise RuntimeError(
                "ConjugationResult does not contain a final Interchange. "
                "Run the build with create_final_interchange enabled before exporting."
            )
        return self.final_interchange

    def get_component_info(self) -> Any:
        """Return component metadata from the retained system builder.

        Returns
        -------
        Any
            Component information returned by ``SystemBuilder.get_component_info()``.

        Raises
        ------
        RuntimeError
            If the result does not retain a compatible system builder.
        """
        if self.system_builder is None or not hasattr(self.system_builder, "get_component_info"):
            raise RuntimeError(
                "ConjugationResult does not contain a system builder with component info. "
                "Component metadata is only available immediately after a full build."
            )
        return self.system_builder.get_component_info()


ConjugateBuildResult = ConjugationResult


def _optional_path(value: Any) -> Path | None:
    if value is None:
        return None
    return Path(value)


def _restore_legacy_pdb_fragment_alias(artifact_paths: dict[str, Path]) -> None:
    """Restore the original unindexed PDB-fragment ingestion artifact key."""
    legacy_key = "pdb_fragment_pdb_fragment_ingestion"
    if legacy_key in artifact_paths:
        return
    indexed_key = "pdb_fragment_1_pdb_fragment_ingestion"
    if indexed_key in artifact_paths:
        artifact_paths[legacy_key] = artifact_paths[indexed_key]
