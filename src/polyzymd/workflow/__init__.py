"""Workflow management for HPC clusters and job submission."""

from polyzymd.workflow.daisy_chain import (
    DaisyChainConfig,
    DaisyChainSubmitter,
    SegmentInfo,
    SubmissionResult,
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
    # Daisy-chain submission
    "DaisyChainConfig",
    "DaisyChainSubmitter",
    "SegmentInfo",
    "SubmissionResult",
    "submit_daisy_chain",
]
