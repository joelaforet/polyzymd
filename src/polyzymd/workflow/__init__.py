"""Workflow management for HPC clusters and job submission."""

from polyzymd.utils.replicates import parse_replicate_range, validate_replicate_range
from polyzymd.workflow.analysis_slurm import (
    AnalysisJobManifest,
    AnalysisSlurmResources,
    ConditionTaskSpec,
    ReplicateTaskSpec,
    SubmittedJobGraph,
    TaskStatus,
    build_manifest,
    generate_aggregate_script,
    generate_finalize_script,
    generate_replicate_script,
    read_analysis_status,
    submit_analysis_graph,
    update_task_status,
)
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
    # Analysis HPC orchestration
    "AnalysisSlurmResources",
    "ReplicateTaskSpec",
    "ConditionTaskSpec",
    "AnalysisJobManifest",
    "SubmittedJobGraph",
    "TaskStatus",
    "build_manifest",
    "generate_replicate_script",
    "generate_aggregate_script",
    "generate_finalize_script",
    "submit_analysis_graph",
    "read_analysis_status",
    "update_task_status",
]
