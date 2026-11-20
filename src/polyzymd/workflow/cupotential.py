"""Site settings for OpenMM jobs that run on NVIDIA's cuPotential plugin.

``openmm.backend: openmm_cupotential`` jobs run in their own pixi workspace,
``envs/openmm-cupotential/pixi.toml``, inside an Apptainer image that provides
glibc >= 2.38 and a CUDA 13.4 driver API. See
``docs/source/how_to/run_with_cupotential.md``.
"""

from __future__ import annotations

import os
from pathlib import Path

#: The pixi environment a cuPotential job records (``runtime_platform.json``).
CUPOTENTIAL_PIXI_ENV = "openmm-cupotential"
IMAGE_ENV_VAR = "POLYZYMD_CUPOTENTIAL_IMAGE"
CHANNEL_ENV_VAR = "POLYZYMD_CUPOTENTIAL_CHANNEL"

#: Blanca node features of GPUs that run the plugin. CUDA 13 dropped V100 and
#: P100, so those nodes are never requested. ``rtx6000`` covers the Quadro RTX
#: 6000 and the RTX PRO 6000.
#: Routing retries for cuPotential jobs. The GPU constraint selects models, not drivers,
#: and on 2026-11-19 about 8 Blanca GPU nodes of those models ran R525/R550 drivers,
#: where cupo-run exits 86. Each bounce excludes one node, so allow enough retries
#: to step past all of them instead of the stock limit of 3.
CUPOTENTIAL_ROUTING_RETRY_LIMIT = 10
STOCK_ROUTING_RETRY_LIMIT = 3

BLANCA_CUPOTENTIAL_FEATURES = ("A40", "A100", "L40", "h100", "rtx6000", "T4")
PRESET_CUPOTENTIAL_CONSTRAINT = {
    "blanca-shirts": "|".join(BLANCA_CUPOTENTIAL_FEATURES),
    "blanca-chbe-rdi": "|".join(BLANCA_CUPOTENTIAL_FEATURES),
}

_SHIRTS_LAB = "/pl/active/shirts_archive/LaforetJoe/containers/openmm-cupotential"
PRESET_CUPOTENTIAL_IMAGE = {"blanca-shirts": f"{_SHIRTS_LAB}/openmm-cupotential_0.3.5.sif"}
PRESET_CUPOTENTIAL_CHANNEL = {
    "blanca-shirts": f"{_SHIRTS_LAB}/conda-channel-0.3.5/nv-openmm-plugin-conda-channel-0.3.5"
}


def resolve_cupotential_image(preset: str | None) -> str:
    """Return the Apptainer image: ``$POLYZYMD_CUPOTENTIAL_IMAGE``, else the preset's.

    Raises
    ------
    ValueError
        When neither names an existing file.
    """
    image = os.environ.get(IMAGE_ENV_VAR) or PRESET_CUPOTENTIAL_IMAGE.get(preset or "")
    if not image:
        raise ValueError(
            f"openmm.backend 'openmm_cupotential' runs in an Apptainer image; set {IMAGE_ENV_VAR} "
            "to the openmm-cupotential .sif file (see docs/source/how_to/run_with_cupotential.md)"
        )
    if not Path(image).is_file():
        raise ValueError(
            f"the openmm-cupotential image {image} does not exist (set {IMAGE_ENV_VAR})"
        )
    return image


def link_cupotential_channel(env_dir: Path, preset: str | None) -> Path:
    """Point ``<env_dir>/nv-channel`` at the plugin's local conda channel and return it.

    The channel is ``$POLYZYMD_CUPOTENTIAL_CHANNEL``, else an existing
    ``nv-channel``, else the preset's. A pixi manifest cannot read environment
    variables, so the manifest names the gitignored link instead.

    Raises
    ------
    ValueError
        When no channel holding ``linux-64/repodata.json`` is found.
    """
    link = env_dir / "nv-channel"
    channel = os.environ.get(CHANNEL_ENV_VAR)
    if not channel and not link.exists():
        channel = PRESET_CUPOTENTIAL_CHANNEL.get(preset or "")
    target = Path(channel) if channel else link
    if not (target / "linux-64" / "repodata.json").is_file():
        raise ValueError(
            f"no nv-openmm-plugin conda channel at {target}: set {CHANNEL_ENV_VAR} to the unpacked "
            "nv-openmm-plugin-conda-channel folder (the one holding linux-64/repodata.json)"
        )
    if channel and link.resolve() != target.resolve():
        if link.is_symlink():
            link.unlink()
        link.symlink_to(target.resolve())
    return link
