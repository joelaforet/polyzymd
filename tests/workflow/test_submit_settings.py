"""Tests for the recorded SLURM submit settings and the site overrides."""

from __future__ import annotations

import pytest

from polyzymd.workflow.slurm import SlurmConfig, sbatch_lines
from polyzymd.workflow.submit_settings import (
    SUBMIT_SETTINGS_FILENAME,
    SubmitSettings,
    apply_slurm_overrides,
    read_submit_settings,
    write_submit_settings,
)


def test_an_empty_qos_or_account_removes_the_preset_line():
    config = apply_slurm_overrides(SlurmConfig.from_preset("aa100"), {"qos": "", "account": ""})

    lines = sbatch_lines(config)
    assert lines["qos_line"] == ""
    assert lines["account_line"] == ""


def test_none_keeps_the_preset_value_and_a_value_replaces_it():
    config = apply_slurm_overrides(
        SlurmConfig.from_preset("aa100"), {"qos": None, "account": "smithlab"}
    )

    assert config.qos == "normal"
    assert config.account == "smithlab"


def test_an_empty_partition_is_refused():
    with pytest.raises(ValueError, match="--partition cannot be empty"):
        apply_slurm_overrides(SlurmConfig.from_preset("aa100"), {"partition": ""})


def test_generic_gpu_preset_sets_no_account_and_no_qos():
    lines = sbatch_lines(SlurmConfig.from_preset("generic-gpu"))

    assert lines["account_line"] == ""
    assert lines["qos_line"] == ""
    assert lines["gpu_line"] == "#SBATCH --gres=gpu:1"


def test_settings_round_trip_through_the_run_folder(tmp_path):
    settings = SubmitSettings.from_options(
        "aa100", pixi_env="sim-cuda-12-0", account="smithlab", qos="", gpu_type="l40s"
    )

    path = write_submit_settings(tmp_path, settings)
    back = read_submit_settings(tmp_path)

    assert path.name == SUBMIT_SETTINGS_FILENAME
    assert back == settings
    config = back.slurm_config()
    assert (config.account, config.qos, config.gpu_type) == ("smithlab", "", "l40s")


def test_command_flags_replay_every_override():
    settings = SubmitSettings.from_options(
        "aa100", pixi_env="sim-cuda-12-0", partition="gpu", account="smithlab", qos=""
    )

    assert settings.command_flags() == (
        "--preset aa100 --partition gpu --account smithlab --qos '' --pixi-env sim-cuda-12-0"
    )


def test_command_flags_add_nodes_to_the_excluded_list():
    preset_nodes = SlurmConfig.from_preset("blanca-shirts").exclude
    recorded = SubmitSettings.from_options("blanca-shirts")
    overridden = SubmitSettings.from_options("blanca-shirts", exclude="n1")

    assert recorded.command_flags(extra_exclude="n9") == (
        f"--preset blanca-shirts --exclude {preset_nodes},n9"
    )
    assert overridden.command_flags(extra_exclude="n9") == "--preset blanca-shirts --exclude n1,n9"


def test_given_options_replace_single_recorded_values():
    recorded = SubmitSettings.from_options("aa100", account="smithlab", qos="")

    merged = recorded.with_options(account="otherlab", time_limit="4:00:00")

    assert merged.overrides == {"account": "otherlab", "qos": "", "time_limit": "4:00:00"}


@pytest.mark.parametrize("text", ["{not json", '{"preset": "nowhere"}', '{"overrides": {}}'])
def test_an_unusable_record_is_ignored(tmp_path, text):
    (tmp_path / SUBMIT_SETTINGS_FILENAME).write_text(text)

    assert read_submit_settings(tmp_path) is None


def test_no_record_gives_none(tmp_path):
    assert read_submit_settings(tmp_path) is None
    assert read_submit_settings(None) is None
