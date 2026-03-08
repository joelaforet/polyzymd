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
    # Job submission
    "DaisyChainConfig",
    "DaisyChainSubmitter",
    "SegmentInfo",  # Deprecated — kept for backward compatibility
    "SubmissionResult",
    "submit_daisy_chain",
]
