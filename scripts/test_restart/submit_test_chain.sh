#!/usr/bin/env bash
# =============================================================================
# Smart-Restart SLURM Daisy Chain Test (Blanca / blanca-shirts)
# =============================================================================
#
# Submits a 3-segment daisy chain on Blanca to test the full smart-restart
# infrastructure under real SLURM conditions.
#
# Chain structure:
#   Segment 0 (setup + initial production) → no dependency
#   Segment 1 (continuation)               → afterany:seg0
#   Segment 2 (continuation + recovery)    → afterany:seg1
#
# Each segment has a 5-minute wall time. The tiny waterbox finishes in seconds,
# so the jobs should complete quickly. The point is to test the SLURM dependency
# chain, signal handling, and recovery preamble.
#
# Interrupt testing:
#   After segment 0 finishes and segment 1 starts, you can manually interrupt
#   segment 1 to test recovery in segment 2:
#
#     scancel --signal=USR1 <segment_1_job_id>
#
#   Segment 2 will detect the INTERRUPTED marker and run `polyzymd recover`
#   before continuing.
#
# Prerequisites:
#   1. On a Blanca login node or with access to blanca-shirts partition
#   2. polyzymd-env conda environment exists
#   3. polyzymd installed from feature/smart-restart branch:
#        mamba run -n polyzymd-env pip install -e ".[dev]"
#   4. This script and setup_test_env.py are in the same directory
#
# Usage:
#   cd /path/to/scripts/test_restart/   # or wherever you want to run
#   bash submit_test_chain.sh
#
# Output:
#   test_restart_workdir/        — simulation data (created by segment 0)
#   slurm_test_seg0.sh           — generated SLURM script for segment 0
#   slurm_test_seg1.sh           — generated SLURM script for segment 1
#   slurm_test_seg2.sh           — generated SLURM script for segment 2
#   slurm-*.out                  — SLURM output logs
#
# Cleanup:
#   rm -rf test_restart_workdir slurm_test_seg*.sh slurm-*.out
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PARTITION="blanca,blanca-shirts"
ACCOUNT="blanca-shirts"
QOS="preemptable"
EXCLUDE="bgpu-bortz1"
TIME_LIMIT="00:05:00"         # 5 minutes per segment
CONDA_ENV="polyzymd-env"
SEGMENT_TIME="0.001"          # 0.001 ns = 500 steps at 2fs timestep
NUM_SAMPLES="10"              # 10 frames per segment

# Resolve paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$(pwd)/test_restart_workdir"
SETUP_SCRIPT="${SCRIPT_DIR}/setup_test_env.py"

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo ""
echo "====================================================================="
echo "  Smart-Restart SLURM Daisy Chain Submission"
echo "====================================================================="
echo ""

if ! command -v sbatch &>/dev/null; then
    echo "ERROR: sbatch not found. Are you on a SLURM cluster?"
    exit 1
fi

if [ ! -f "$SETUP_SCRIPT" ]; then
    echo "ERROR: setup_test_env.py not found at $SETUP_SCRIPT"
    echo "       Make sure you're running from the scripts/test_restart/ directory"
    echo "       or that setup_test_env.py is in the same directory as this script."
    exit 1
fi

echo "Configuration:"
echo "  Partition:    $PARTITION"
echo "  Account:      $ACCOUNT"
echo "  QoS:          $QOS"
echo "  Exclude:      $EXCLUDE"
echo "  Wall time:    $TIME_LIMIT (per segment)"
echo "  Conda env:    $CONDA_ENV"
echo "  Segment time: $SEGMENT_TIME ns"
echo "  Num samples:  $NUM_SAMPLES"
echo "  Working dir:  $WORKDIR"
echo ""

# ---------------------------------------------------------------------------
# Generate SLURM scripts
# ---------------------------------------------------------------------------

# --- Segment 0: Setup + initial production ---------------------------------
cat > slurm_test_seg0.sh <<'SLURM_SEG0'
#!/bin/bash
#SBATCH --partition=PARTITION_PLACEHOLDER
#SBATCH --account=ACCOUNT_PLACEHOLDER
#SBATCH --qos=QOS_PLACEHOLDER
#SBATCH --exclude=EXCLUDE_PLACEHOLDER
#SBATCH --job-name=test_seg0
#SBATCH --output=slurm-seg0-%j.out
#SBATCH --time=TIME_PLACEHOLDER
#SBATCH --ntasks=1
#SBATCH --signal=B:USR1@120
#SBATCH --no-requeue

# =============================================================================
# Segment 0: Setup waterbox + run initial production
# =============================================================================

module purge 2>/dev/null || true
ml miniforge 2>/dev/null || true

eval "$(conda shell.bash hook)"
mamba activate CONDA_ENV_PLACEHOLDER

set -e

echo "===== Segment 0: Setup + Initial Production ====="
echo "Hostname: $(hostname)"
echo "Timestamp: $(date)"
echo ""

# Run setup script (creates workdir with waterbox + production_0)
cd SCRIPT_DIR_PLACEHOLDER
python setup_test_env.py

echo ""
echo "Segment 0 complete at $(date)"
SLURM_SEG0

# Replace placeholders
sed -i "s|PARTITION_PLACEHOLDER|${PARTITION}|g" slurm_test_seg0.sh
sed -i "s|ACCOUNT_PLACEHOLDER|${ACCOUNT}|g" slurm_test_seg0.sh
sed -i "s|QOS_PLACEHOLDER|${QOS}|g" slurm_test_seg0.sh
sed -i "s|EXCLUDE_PLACEHOLDER|${EXCLUDE}|g" slurm_test_seg0.sh
sed -i "s|TIME_PLACEHOLDER|${TIME_LIMIT}|g" slurm_test_seg0.sh
sed -i "s|CONDA_ENV_PLACEHOLDER|${CONDA_ENV}|g" slurm_test_seg0.sh
sed -i "s|SCRIPT_DIR_PLACEHOLDER|${SCRIPT_DIR}|g" slurm_test_seg0.sh

# --- Segment 1: Continuation (user can scancel --signal=USR1 this one) ------
cat > slurm_test_seg1.sh <<'SLURM_SEG1'
#!/bin/bash
#SBATCH --partition=PARTITION_PLACEHOLDER
#SBATCH --account=ACCOUNT_PLACEHOLDER
#SBATCH --qos=QOS_PLACEHOLDER
#SBATCH --exclude=EXCLUDE_PLACEHOLDER
#SBATCH --job-name=test_seg1
#SBATCH --output=slurm-seg1-%j.out
#SBATCH --time=TIME_PLACEHOLDER
#SBATCH --ntasks=1
#SBATCH --signal=B:USR1@120
#SBATCH --no-requeue

# =============================================================================
# Segment 1: Continuation from segment 0
# This is the segment to interrupt for testing:
#   scancel --signal=USR1 <this_job_id>
# =============================================================================

module purge 2>/dev/null || true
ml miniforge 2>/dev/null || true

eval "$(conda shell.bash hook)"
mamba activate CONDA_ENV_PLACEHOLDER

echo "===== Segment 1: Continuation from Segment 0 ====="
echo "Hostname: $(hostname)"
echo "Timestamp: $(date)"
echo ""

# ---------------------------------------------------
# Recovery preamble: check previous segment status
# ---------------------------------------------------
PREV_SEG=0
PREV_DIR="WORKDIR_PLACEHOLDER/production_${PREV_SEG}"

if [ -f "${PREV_DIR}/INTERRUPTED" ]; then
    echo "Previous segment $PREV_SEG was interrupted — running recovery"
    polyzymd recover -w "WORKDIR_PLACEHOLDER" -s "$PREV_SEG"
    RECOVER_RC=$?
    if [ $RECOVER_RC -ne 0 ]; then
        echo "Recovery of segment $PREV_SEG FAILED (exit code $RECOVER_RC)"
        exit 1
    fi
    echo "Recovery of segment $PREV_SEG completed successfully"
elif [ -f "${PREV_DIR}/production_${PREV_SEG}_state.xml" ]; then
    echo "Previous segment $PREV_SEG completed normally — proceeding"
else
    echo "ERROR: Previous segment $PREV_SEG has no state.xml and no INTERRUPTED marker"
    echo "This segment cannot continue — the previous job likely crashed."
    exit 1
fi

set -e

echo ""
echo "Running continuation for segment 1..."

cd SCRIPT_DIR_PLACEHOLDER

polyzymd continue \
    -w "WORKDIR_PLACEHOLDER" \
    -s 1 \
    -t SEGMENT_TIME_PLACEHOLDER \
    -n NUM_SAMPLES_PLACEHOLDER
RC=$?

# Exit code 99 = interrupted but state saved (graceful shutdown)
if [ $RC -eq 99 ]; then
    echo "Segment 1 interrupted (graceful shutdown) at $(date)"
    exit 99
elif [ $RC -ne 0 ]; then
    echo "Segment 1 FAILED with exit code $RC at $(date)"
    exit $RC
fi

echo ""
echo "Segment 1 completed successfully at $(date)"
SLURM_SEG1

sed -i "s|PARTITION_PLACEHOLDER|${PARTITION}|g" slurm_test_seg1.sh
sed -i "s|ACCOUNT_PLACEHOLDER|${ACCOUNT}|g" slurm_test_seg1.sh
sed -i "s|QOS_PLACEHOLDER|${QOS}|g" slurm_test_seg1.sh
sed -i "s|EXCLUDE_PLACEHOLDER|${EXCLUDE}|g" slurm_test_seg1.sh
sed -i "s|TIME_PLACEHOLDER|${TIME_LIMIT}|g" slurm_test_seg1.sh
sed -i "s|CONDA_ENV_PLACEHOLDER|${CONDA_ENV}|g" slurm_test_seg1.sh
sed -i "s|SCRIPT_DIR_PLACEHOLDER|${SCRIPT_DIR}|g" slurm_test_seg1.sh
sed -i "s|WORKDIR_PLACEHOLDER|${WORKDIR}|g" slurm_test_seg1.sh
sed -i "s|SEGMENT_TIME_PLACEHOLDER|${SEGMENT_TIME}|g" slurm_test_seg1.sh
sed -i "s|NUM_SAMPLES_PLACEHOLDER|${NUM_SAMPLES}|g" slurm_test_seg1.sh

# --- Segment 2: Continuation with recovery (auto-recovers from segment 1) ---
cat > slurm_test_seg2.sh <<'SLURM_SEG2'
#!/bin/bash
#SBATCH --partition=PARTITION_PLACEHOLDER
#SBATCH --account=ACCOUNT_PLACEHOLDER
#SBATCH --qos=QOS_PLACEHOLDER
#SBATCH --exclude=EXCLUDE_PLACEHOLDER
#SBATCH --job-name=test_seg2
#SBATCH --output=slurm-seg2-%j.out
#SBATCH --time=TIME_PLACEHOLDER
#SBATCH --ntasks=1
#SBATCH --signal=B:USR1@120
#SBATCH --no-requeue

# =============================================================================
# Segment 2: Continuation with recovery preamble
# If segment 1 was interrupted (via scancel --signal=USR1), this segment
# will automatically detect the INTERRUPTED marker, run `polyzymd recover`
# to finish segment 1's remaining steps, then continue with segment 2.
# =============================================================================

module purge 2>/dev/null || true
ml miniforge 2>/dev/null || true

eval "$(conda shell.bash hook)"
mamba activate CONDA_ENV_PLACEHOLDER

echo "===== Segment 2: Continuation with Recovery Preamble ====="
echo "Hostname: $(hostname)"
echo "Timestamp: $(date)"
echo ""

# ---------------------------------------------------
# Recovery preamble: check previous segment status
# ---------------------------------------------------
PREV_SEG=1
PREV_DIR="WORKDIR_PLACEHOLDER/production_${PREV_SEG}"

if [ -f "${PREV_DIR}/INTERRUPTED" ]; then
    echo "Previous segment $PREV_SEG was interrupted — running recovery"
    echo "INTERRUPTED marker contents:"
    cat "${PREV_DIR}/INTERRUPTED"
    echo ""

    polyzymd recover -w "WORKDIR_PLACEHOLDER" -s "$PREV_SEG"
    RECOVER_RC=$?
    if [ $RECOVER_RC -ne 0 ]; then
        echo "Recovery of segment $PREV_SEG FAILED (exit code $RECOVER_RC)"
        exit 1
    fi
    echo "Recovery of segment $PREV_SEG completed successfully"

    # Verify recovery produced the required state file
    if [ ! -f "${PREV_DIR}/production_${PREV_SEG}_state.xml" ]; then
        echo "ERROR: Recovery did not produce state.xml for segment $PREV_SEG"
        exit 1
    fi
    echo "Verified: production_${PREV_SEG}_state.xml exists after recovery"
elif [ -f "${PREV_DIR}/production_${PREV_SEG}_state.xml" ]; then
    echo "Previous segment $PREV_SEG completed normally — proceeding"
else
    echo "ERROR: Previous segment $PREV_SEG has no state.xml and no INTERRUPTED marker"
    echo "This segment cannot continue — the previous job likely crashed."
    exit 1
fi

set -e

echo ""
echo "Running continuation for segment 2..."

cd SCRIPT_DIR_PLACEHOLDER

polyzymd continue \
    -w "WORKDIR_PLACEHOLDER" \
    -s 2 \
    -t SEGMENT_TIME_PLACEHOLDER \
    -n NUM_SAMPLES_PLACEHOLDER
RC=$?

if [ $RC -eq 99 ]; then
    echo "Segment 2 interrupted (graceful shutdown) at $(date)"
    exit 99
elif [ $RC -ne 0 ]; then
    echo "Segment 2 FAILED with exit code $RC at $(date)"
    exit $RC
fi

echo ""
echo "Segment 2 completed successfully at $(date)"

# ---------------------------------------------------
# Final status report
# ---------------------------------------------------
echo ""
echo "===== Chain Complete: Final Directory Status ====="
for d in WORKDIR_PLACEHOLDER/production_*; do
    if [ -d "$d" ]; then
        name=$(basename "$d")
        n_files=$(ls "$d" | wc -l)
        status="UNKNOWN"
        if [ -f "$d/INTERRUPTED" ]; then
            status="INTERRUPTED"
        elif ls "$d"/*_state.xml &>/dev/null 2>&1; then
            status="COMPLETED"
        else
            status="MISSING STATE"
        fi
        echo "  $name/ — $n_files files — $status"
    fi
done
echo ""
SLURM_SEG2

sed -i "s|PARTITION_PLACEHOLDER|${PARTITION}|g" slurm_test_seg2.sh
sed -i "s|ACCOUNT_PLACEHOLDER|${ACCOUNT}|g" slurm_test_seg2.sh
sed -i "s|QOS_PLACEHOLDER|${QOS}|g" slurm_test_seg2.sh
sed -i "s|EXCLUDE_PLACEHOLDER|${EXCLUDE}|g" slurm_test_seg2.sh
sed -i "s|TIME_PLACEHOLDER|${TIME_LIMIT}|g" slurm_test_seg2.sh
sed -i "s|CONDA_ENV_PLACEHOLDER|${CONDA_ENV}|g" slurm_test_seg2.sh
sed -i "s|SCRIPT_DIR_PLACEHOLDER|${SCRIPT_DIR}|g" slurm_test_seg2.sh
sed -i "s|WORKDIR_PLACEHOLDER|${WORKDIR}|g" slurm_test_seg2.sh
sed -i "s|SEGMENT_TIME_PLACEHOLDER|${SEGMENT_TIME}|g" slurm_test_seg2.sh
sed -i "s|NUM_SAMPLES_PLACEHOLDER|${NUM_SAMPLES}|g" slurm_test_seg2.sh

echo "Generated SLURM scripts:"
echo "  slurm_test_seg0.sh"
echo "  slurm_test_seg1.sh"
echo "  slurm_test_seg2.sh"
echo ""

# ---------------------------------------------------------------------------
# Submit the daisy chain
# ---------------------------------------------------------------------------
echo "Submitting segment 0 (setup + initial production)..."
SEG0_RESULT=$(sbatch slurm_test_seg0.sh)
SEG0_ID=$(echo "$SEG0_RESULT" | awk '{print $NF}')
echo "  $SEG0_RESULT"
echo ""

echo "Submitting segment 1 (continuation, afterany:$SEG0_ID)..."
SEG1_RESULT=$(sbatch --dependency=afterany:${SEG0_ID} slurm_test_seg1.sh)
SEG1_ID=$(echo "$SEG1_RESULT" | awk '{print $NF}')
echo "  $SEG1_RESULT"
echo ""

echo "Submitting segment 2 (continuation + recovery, afterany:$SEG1_ID)..."
SEG2_RESULT=$(sbatch --dependency=afterany:${SEG1_ID} slurm_test_seg2.sh)
SEG2_ID=$(echo "$SEG2_RESULT" | awk '{print $NF}')
echo "  $SEG2_RESULT"
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "====================================================================="
echo "  Daisy Chain Submitted Successfully"
echo "====================================================================="
echo ""
echo "  Segment 0 (setup):        Job $SEG0_ID"
echo "  Segment 1 (continue):     Job $SEG1_ID  (afterany:$SEG0_ID)"
echo "  Segment 2 (continue+rec): Job $SEG2_ID  (afterany:$SEG1_ID)"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  watch -n5 squeue -u \$USER"
echo ""
echo "To test signal-based interruption, wait for segment 1 to start, then:"
echo "  scancel --signal=USR1 $SEG1_ID"
echo ""
echo "This will:"
echo "  1. Send SIGUSR1 to segment 1's Python process"
echo "  2. Trigger graceful shutdown (save emergency state, exit 99)"
echo "  3. Segment 2 will detect the INTERRUPTED marker and auto-recover"
echo ""
echo "Check results after all jobs complete:"
echo "  cat slurm-seg0-${SEG0_ID}.out"
echo "  cat slurm-seg1-${SEG1_ID}.out"
echo "  cat slurm-seg2-${SEG2_ID}.out"
echo "  ls -la test_restart_workdir/production_*/"
echo ""
echo "If segment 1 was NOT interrupted, all 3 segments complete normally."
echo "If segment 1 WAS interrupted, segment 2 recovers it automatically."
echo ""
echo "Cleanup:"
echo "  rm -rf test_restart_workdir slurm_test_seg*.sh slurm-seg*-*.out"
echo ""
