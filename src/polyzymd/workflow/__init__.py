"""Workflow management for HPC clusters and job submission."""

from polyzymd.workflow.daisy_chain import (
    DaisyChainConfig,
    DaisyChainSubmitter,
    SubmissionResult,
    create_job_name,
    submit_daisy_chain,
)
from polyzymd.workflow.slurm import (
    JobContext,
    SlurmConfig,
    SlurmScriptGenerator,
    parse_replicate_range,
    validate_replicate_range,
)

__all__ = [
    # SLURM utilities
    "JobContext",
    "SlurmConfig",
    "SlurmScriptGenerator",
    "parse_replicate_range",
    "validate_replicate_range",
    # Job submission
    "DaisyChainConfig",
    "DaisyChainSubmitter",
    "SubmissionResult",
    "create_job_name",
    "submit_daisy_chain",
]
