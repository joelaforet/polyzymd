"""Conjugation orchestration boundary."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from polyzymd.builders.conjugation.models import ConjugateBuildRequest, ConjugationResult

if TYPE_CHECKING:
    from polyzymd.builders.conjugation.system_workflow import ConjugatedPolymerSystemSettings


class ConjugationEngine:
    """Boundary for public conjugate construction orchestration.

    The engine currently delegates to the existing config-driven conjugated
    polymer workflow. Future phases should move orchestration behind this class
    without duplicating behavior in facade functions.
    """

    def __init__(self, *, settings: ConjugatedPolymerSystemSettings | None = None) -> None:
        self.settings = settings

    def build_from_config(
        self,
        config: Any | Path | str,
        *,
        output_dir: Path | str | None = None,
        settings: ConjugatedPolymerSystemSettings | None = None,
        free_polymer_seed: int | None = None,
    ) -> ConjugationResult:
        """Delegate config-based construction to the existing workflow."""
        from polyzymd.builders.conjugation.system_workflow import (
            build_conjugated_polymer_system_from_config,
            build_conjugated_polymer_system_from_config_path,
        )

        workflow_settings = self.settings if settings is None else settings
        if isinstance(config, (str, Path)):
            legacy_result = build_conjugated_polymer_system_from_config_path(
                config,
                output_dir=output_dir,
                settings=workflow_settings,
                free_polymer_seed=free_polymer_seed,
            )
            return ConjugationResult.from_legacy_result(legacy_result, config_path=config)

        if output_dir is None:
            raise ValueError(
                "output_dir is required when build_from_config() receives an "
                "in-memory SimulationConfig."
            )
        legacy_result = build_conjugated_polymer_system_from_config(
            config,
            output_dir=output_dir,
            settings=workflow_settings,
            free_polymer_seed=free_polymer_seed,
        )
        return ConjugationResult.from_legacy_result(legacy_result)

    def build_from_request(
        self,
        request: ConjugateBuildRequest,
        *,
        settings: ConjugatedPolymerSystemSettings | None = None,
    ) -> ConjugationResult:
        """Build from the lightweight public request shell."""
        source = request.config if request.config is not None else request.config_path
        if source is None:
            return self._build_direct_request(request, settings=settings)
        return self.build_from_config(
            source,
            output_dir=request.output_dir,
            settings=settings,
            free_polymer_seed=request.free_polymer_seed,
        )

    def build(
        self,
        request: ConjugateBuildRequest | Any | Path | str | None = None,
        **kwargs: Any,
    ) -> ConjugationResult:
        """Build from a request, config object, or config path.

        Direct molecule/topology construction is intentionally not routed to a
        silent OpenFF fallback. Use the config-driven path until a future engine
        phase owns that workflow explicitly.
        """
        call_settings = kwargs.pop("settings", None)

        if isinstance(request, ConjugateBuildRequest):
            if kwargs:
                _raise_unsupported_direct_mode(kwargs)
            return self.build_from_request(request, settings=call_settings)

        if request is None:
            if "config" in kwargs:
                request = kwargs.pop("config")
            elif "config_path" in kwargs:
                request = kwargs.pop("config_path")
            elif "protein_pdb_path" in kwargs or "attachments" in kwargs:
                request = ConjugateBuildRequest(**kwargs)
                return self.build_from_request(request, settings=call_settings)
            else:
                _raise_unsupported_direct_mode(kwargs)

        output_dir = kwargs.pop("output_dir", None)
        free_polymer_seed = kwargs.pop("free_polymer_seed", None)
        if kwargs:
            _raise_unsupported_direct_mode(kwargs)

        return self.build_from_config(
            request,
            output_dir=output_dir,
            settings=call_settings,
            free_polymer_seed=free_polymer_seed,
        )

    def _build_direct_request(
        self,
        request: ConjugateBuildRequest,
        *,
        settings: ConjugatedPolymerSystemSettings | None = None,
    ) -> ConjugationResult:
        """Build a direct protein plus attachment request."""
        if request.protein_pdb_path is None or not request.attachments:
            raise ValueError(
                "ConjugateBuildRequest requires config/config_path or protein_pdb_path with attachments"
            )
        if request.output_dir is None:
            raise ValueError("output_dir is required for direct conjugation requests")

        from polyzymd.builders.conjugation.system_workflow import (
            build_direct_smiles_moiety_conjugate,
        )

        workflow_settings = self.settings if settings is None else settings
        legacy_result = build_direct_smiles_moiety_conjugate(
            protein_pdb_path=request.protein_pdb_path,
            attachments=request.attachments,
            output_dir=request.output_dir,
            ccd_pablo=request.ccd_pablo,
            chain_policy=request.chain_policy,
            settings=workflow_settings,
            random_seed=request.free_polymer_seed,
        )
        return ConjugationResult.from_legacy_result(legacy_result)


def _raise_unsupported_direct_mode(inputs: dict[str, Any]) -> None:
    supplied = ", ".join(sorted(inputs)) if inputs else "no config input"
    raise NotImplementedError(
        "Direct ConjugationEngine.build(...) molecule/topology inputs are not implemented. "
        "Use config, config_path, or ConjugateBuildRequest to delegate to the existing "
        f"config-driven workflow. Received: {supplied}."
    )
