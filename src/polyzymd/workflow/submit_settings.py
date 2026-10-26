"""The SLURM settings a replicate was submitted with.

``polyzymd submit`` takes a preset and site overrides (account, QoS,
partition, GPU type, ...). It records them in each replicate's run folder
(:data:`SUBMIT_SETTINGS_FILENAME`), so that ``polyzymd recover --submit`` and
the resubmit commands of ``polyzymd status --format agent`` reuse them instead
of the preset's values. A replicate without the record keeps the preset-only
behaviour.
"""

from __future__ import annotations

import json
import logging
import shlex
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from polyzymd.workflow.slurm import SlurmConfig

LOGGER = logging.getLogger(__name__)

SUBMIT_SETTINGS_FILENAME = "submit_settings.json"

#: SlurmConfig fields that ``submit`` and ``recover`` can override, in flag order.
OVERRIDE_FIELDS = (
    "partition",
    "account",
    "qos",
    "time_limit",
    "memory",
    "gpu_type",
    "constraint",
    "nodelist",
    "exclude",
    "email",
)
# Fields whose empty value removes the ``#SBATCH`` line. A job always needs a
# partition and a time limit, so those two cannot be emptied.
_REMOVABLE = {
    "account": "",
    "qos": "",
    "email": "",
    "memory": None,
    "gpu_type": None,
    "constraint": None,
    "nodelist": None,
    "exclude": None,
}


def apply_slurm_overrides(config: SlurmConfig, overrides: Mapping[str, str | None]) -> SlurmConfig:
    """Apply command-line overrides to ``config`` in place and return it.

    ``None`` keeps the preset's value. An empty string removes the setting,
    so its ``#SBATCH`` line is not written (``--qos ""`` on a site without
    QoS).

    Raises
    ------
    ValueError
        For an unknown field, or an empty partition or time limit.
    """
    for key, value in overrides.items():
        if key not in OVERRIDE_FIELDS:
            raise ValueError(f"unknown SLURM override {key!r}")
        if value is None:
            continue
        if value == "":
            if key not in _REMOVABLE:
                raise ValueError(f"--{key.replace('_', '-')} cannot be empty")
            setattr(config, key, _REMOVABLE[key])
        else:
            setattr(config, key, value)
    return config


@dataclass
class SubmitSettings:
    """A preset, the overrides given with it, and the job's pixi environment.

    ``overrides`` holds only the settings given on the command line; an
    empty string is a removed setting.
    """

    preset: str
    overrides: dict[str, str] = field(default_factory=dict)
    pixi_env: str | None = None

    @classmethod
    def from_options(
        cls, preset: str, pixi_env: str | None = None, **options: str | None
    ) -> SubmitSettings:
        """Keep the options that were given (not ``None``)."""
        overrides = {key: value for key, value in options.items() if value is not None}
        apply_slurm_overrides(SlurmConfig(), overrides)  # refuse unknown or empty fields now
        return cls(preset=preset, overrides=overrides, pixi_env=pixi_env)

    def slurm_config(self) -> SlurmConfig:
        """The preset with the overrides applied."""
        return apply_slurm_overrides(SlurmConfig.from_preset(self.preset), self.overrides)  # type: ignore[arg-type]

    def with_options(self, pixi_env: str | None = None, **options: str | None) -> SubmitSettings:
        """These settings with the given options replacing the recorded ones."""
        overrides = dict(self.overrides)
        overrides.update({key: value for key, value in options.items() if value is not None})
        return SubmitSettings.from_options(
            self.preset, pixi_env=pixi_env or self.pixi_env, **overrides
        )

    def command_flags(self, extra_exclude: str = "") -> str:
        """The ``submit``/``recover`` flags that give these settings, shell-quoted.

        ``extra_exclude`` adds nodes to the excluded list, since ``--exclude``
        replaces the preset's list rather than adding to it.
        """
        overrides = dict(self.overrides)
        if extra_exclude:
            base = overrides.get("exclude")
            if base is None:
                base = SlurmConfig.from_preset(self.preset).exclude or ""  # type: ignore[arg-type]
            names = [n for n in f"{base},{extra_exclude}".split(",") if n]
            overrides["exclude"] = ",".join(dict.fromkeys(names))
        words = ["--preset", self.preset]
        for key in OVERRIDE_FIELDS:
            if key in overrides:
                words += [f"--{key.replace('_', '-')}", overrides[key]]
        if self.pixi_env:
            words += ["--pixi-env", self.pixi_env]
        return " ".join(shlex.quote(word) for word in words)

    def to_dict(self) -> dict:
        return {"preset": self.preset, "overrides": self.overrides, "pixi_env": self.pixi_env}


def write_submit_settings(run_dir: Path | str, settings: SubmitSettings) -> Path:
    """Record ``settings`` in ``run_dir`` and return the file."""
    path = Path(run_dir) / SUBMIT_SETTINGS_FILENAME
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(settings.to_dict(), indent=1) + "\n")
    tmp.replace(path)
    return path


def read_submit_settings(run_dir: Path | str | None) -> SubmitSettings | None:
    """Return the settings recorded in ``run_dir``, or ``None`` without a usable record."""
    if run_dir is None:
        return None
    path = Path(run_dir) / SUBMIT_SETTINGS_FILENAME
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text())
        overrides = {str(k): str(v) for k, v in dict(data.get("overrides") or {}).items()}
        settings = SubmitSettings.from_options(
            str(data["preset"]), pixi_env=data.get("pixi_env"), **overrides
        )
        settings.slurm_config()  # an unknown preset is not a usable record
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        LOGGER.warning(f"Ignoring {path}: {exc}")
        return None
    return settings
