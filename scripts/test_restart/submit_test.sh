#!/usr/bin/env bash
# =============================================================================
# Self-Resubmitting SLURM Job Test (Blanca / blanca-shirts)
# =============================================================================
#
# Submits a SINGLE self-resubmitting job on Blanca that tests the full
# self-resubmitting architecture under real SLURM conditions.
#
# How it works:
#   1. setup_test_env.py creates a tiny waterbox + segment 0 + config + progress
#   2. This script generates a self-resubmitting SLURM script
#   3. The job calls `polyzymd run-segment` which:
#      - Loads config + progress
#      - Determines what work remains
#      - Runs the next segment
#   4. After the segment, it calls `polyzymd check-progress`:
#      - Exit 0 = complete, job ends
#      - Exit 1 = work remains, job resubmits itself via `sbatch "$SLURM_JOB_SCRIPT"`
#
# Signal testing:
#   While the job is running, you can interrupt it with:
#     scancel --signal=USR1 <job_id>
#   The job will save interrupted state (exit 99), then resubmit. The next
#   invocation will detect the interrupted segment and continue from there.
#
# Prerequisites:
#   1. On a Blanca login node or with access to blanca-shirts partition
#   2. A PolyzyMD pixi environment exists
#   3. PolyzyMD is installed via `pixi install -e <env>`
#   4. setup_test_env.py has already been run:
#        pixi run -e cuda-12-4 python scripts/test_restart/setup_test_env.py
#
# Usage:
#   cd /path/to/scripts/test_restart/
#   pixi run -e cuda-12-4 python scripts/test_restart/setup_test_env.py   # if not already done
#   bash submit_test.sh
#
# Output:
#   slurm_test_selfresubmit.sh   — the generated SLURM script
#   slurm-selfresubmit-*.out     — SLURM output logs
#
# Cleanup:
#   rm -rf test_restart_workdir test_config.yaml slurm_test_selfresubmit.sh slurm-selfresubmit-*.out
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PARTITION="blanca,blanca-shirts"
ACCOUNT="blanca-shirts"
QOS="preemptable"
EXCLUDE="bgpu-bortz1"
WALL_TIME="00:05:00"          # 5 minutes per job invocation
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIXI_ENV="${PIXI_ENV:-cuda-12-4}"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PIXIMANIFEST="${REPO_ROOT}/pixi.toml"

# Resolve paths
WORKDIR="${SCRIPT_DIR}/test_restart_workdir"
CONFIG="${SCRIPT_DIR}/test_config.yaml"
SLURM_SCRIPT="${SCRIPT_DIR}/slurm_test_selfresubmit.sh"

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo ""
echo "====================================================================="
echo "  Self-Resubmitting SLURM Job Test"
echo "====================================================================="
echo ""

if ! command -v sbatch &>/dev/null; then
    echo "ERROR: sbatch not found. Are you on a SLURM cluster?"
    exit 1
fi

if [ ! -d "$WORKDIR" ]; then
    echo "ERROR: $WORKDIR not found."
    echo "       Run setup_test_env.py first:"
    echo "         pixi run -e ${PIXI_ENV} python scripts/test_restart/setup_test_env.py"
    exit 1
fi

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: $CONFIG not found."
    echo "       Run setup_test_env.py first."
    exit 1
fi

if [ ! -f "$WORKDIR/progress.json" ]; then
    echo "ERROR: progress.json not found in $WORKDIR"
    echo "       Run setup_test_env.py first."
    exit 1
fi

echo "Configuration:"
echo "  Partition:    $PARTITION"
echo "  Account:      $ACCOUNT"
echo "  QoS:          $QOS"
echo "  Exclude:      $EXCLUDE"
echo "  Wall time:    $WALL_TIME (per job invocation)"
echo "  Pixi env:     $PIXI_ENV"
echo "  Working dir:  $WORKDIR"
echo "  Config:       $CONFIG"
echo ""

# ---------------------------------------------------------------------------
# Generate the self-resubmitting SLURM script
# ---------------------------------------------------------------------------
cat > "$SLURM_SCRIPT" <<SLURM_EOF
#!/bin/bash
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --qos=${QOS}
#SBATCH --exclude=${EXCLUDE}
#SBATCH --job-name=test_selfresubmit
#SBATCH --output=${SCRIPT_DIR}/slurm-selfresubmit-%j.out
#SBATCH --time=${WALL_TIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --gres=gpu:1
#SBATCH --signal=B:USR1@120
#SBATCH --no-requeue

# =============================================================================
# Self-Resubmitting Job: waterbox restart test
#
# This script is submitted once and resubmits itself until the simulation
# is complete. Each invocation calls \`polyzymd run-segment\` which
# automatically determines what work remains.
# =============================================================================

echo "===== Self-Resubmitting Job: Invocation at \$(date) ====="
echo "Job ID: \$SLURM_JOB_ID"
echo "Hostname: \$(hostname)"
echo ""

# --- Environment setup ---
eval "\$(pixi shell-hook -e ${PIXI_ENV} --manifest-path ${PIXIMANIFEST})"

# --- Signal forwarding (trap + background + wait) ---
CHILD_PID=""
forward_signal() {
    if [ -n "\$CHILD_PID" ] && kill -0 "\$CHILD_PID" 2>/dev/null; then
        echo "Forwarding \$1 to Python process (PID \$CHILD_PID)"
        kill -"\$1" "\$CHILD_PID"
    fi
}
trap 'forward_signal USR1' USR1
trap 'forward_signal TERM' TERM

# --- Run the next segment ---
echo "Running: polyzymd run-segment -c ${CONFIG} -r 1 --scratch-dir ${WORKDIR}"
echo ""

polyzymd run-segment \\
    -c "${CONFIG}" \\
    -r 1 \\
    --scratch-dir "${WORKDIR}" &
CHILD_PID=\$!

# Wait for the child; disable set -e to capture non-zero exit codes
set +e
wait "\$CHILD_PID" 2>/dev/null
RC=\$?
while kill -0 "\$CHILD_PID" 2>/dev/null; do
    wait "\$CHILD_PID" 2>/dev/null
    RC=\$?
done
set -e

echo ""
echo "run-segment exited with code: \$RC"

# --- Handle exit codes ---
if [ \$RC -eq 1 ]; then
    echo "ERROR: run-segment failed. Check logs above."
    exit 1
fi

# Exit 0 (completed segment) or 99 (interrupted) — check if more work remains
echo ""
echo "Checking progress..."
set +e
polyzymd check-progress \\
    -c "${CONFIG}" \\
    -r 1 \\
    --scratch-dir "${WORKDIR}"
PROGRESS_RC=\$?
set -e

if [ \$PROGRESS_RC -eq 0 ]; then
    echo ""
    echo "===== Simulation COMPLETE at \$(date) ====="
    echo "No resubmission needed."
    exit 0
fi

# --- Resubmit ---
# NOTE: \$SLURM_JOB_SCRIPT is only available in SLURM >= 22.05.
# Use the hardcoded script path as the reliable fallback.
SCRIPT_PATH="${SLURM_SCRIPT}"
echo ""
echo "Work remains. Resubmitting..."
echo "  Script: \$SCRIPT_PATH"
RESUB_RESULT=\$(sbatch "\$SCRIPT_PATH")
echo "Resubmitted: \$RESUB_RESULT"
echo ""
echo "===== Job ending at \$(date), continuation submitted ====="
SLURM_EOF

echo "Generated SLURM script: $SLURM_SCRIPT"
echo ""

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
echo "Submitting self-resubmitting job..."
JOB_RESULT=$(sbatch "$SLURM_SCRIPT")
JOB_ID=$(echo "$JOB_RESULT" | awk '{print $NF}')
echo "  $JOB_RESULT"
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "====================================================================="
echo "  Job Submitted Successfully"
echo "====================================================================="
echo ""
echo "  Job ID: $JOB_ID"
echo "  Script: $SLURM_SCRIPT"
echo ""
echo "  The job will run segments until the simulation is complete,"
echo "  resubmitting itself after each segment."
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  watch -n5 squeue -u \$USER"
echo ""
echo "To test signal-based interruption (while job is running):"
echo "  scancel --batch --signal=USR1 $JOB_ID"
echo ""
echo "  NOTE: --batch is required to send the signal to the batch shell"
echo "  (where our trap + forward_signal logic runs). Without --batch,"
echo "  scancel sends to job steps, which bypasses the batch script."
echo ""
echo "This will:"
echo "  1. Send SIGUSR1 to the Python process"
echo "  2. Trigger graceful shutdown (save interrupted state, exit 99)"
echo "  3. Job resubmits itself"
echo "  4. Next invocation auto-recovers from the interrupted segment"
echo ""
echo "View logs as they appear:"
echo "  tail -f ${SCRIPT_DIR}/slurm-selfresubmit-*.out"
echo ""
echo "Check progress at any time:"
echo "  pixi run -e $PIXI_ENV polyzymd check-progress -c $CONFIG -r 1 --scratch-dir $WORKDIR"
echo ""
echo "Cleanup when done:"
echo "  rm -rf $WORKDIR $CONFIG $SLURM_SCRIPT ${SCRIPT_DIR}/slurm-selfresubmit-*.out"
echo ""
