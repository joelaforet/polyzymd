"""Engine-neutral enhanced sampling placeholders."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


class ExternalBiasSpec(BaseModel):
    """External bias specification.

    Planned for future release. Not yet implemented.

    Parameters
    ----------
    bias_type : str, optional
        Bias backend name, by default ``"plumed"``.
    script_path : Path | None, optional
        Path to an external bias script.
    script_content : str | None, optional
        Inline script content when no file path is used.
    """

    bias_type: str = "plumed"
    script_path: Path | None = None
    script_content: str | None = None


class EnhancedSamplingProtocol(BaseModel):
    """Enhanced sampling protocol definition.

    Planned for future release. Not yet implemented.

    Parameters
    ----------
    protocol_type : str
        Enhanced-sampling protocol identifier.
    parameters : dict[str, Any], optional
        Protocol-specific parameter mapping.
    """

    protocol_type: str
    parameters: dict[str, Any] = Field(default_factory=dict)
