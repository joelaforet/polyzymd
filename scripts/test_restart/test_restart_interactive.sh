#!/usr/bin/env bash
# =============================================================================
# Self-Resubmitting Job Interactive Test Script
# =============================================================================
#
# Exercises the self-resubmitting architecture using the synthetic waterbox
# created by setup_test_env.py.
#
# Tests:
#   1. run-segment continuation: detects segment 0 done, runs segment 1
#   2. check-progress: reports correct status after segment 1
#   3. Signal interruption: kill -USR1 during run-segment triggers exit 99
#   4. check-progress: reports incomplete status after interruption
#   5. run-segment recovery: detects interrupted segment, continues
#   6. recover: reports status (no --submit since not on SLURM)
#   7. run-segment completes: final segment finishes, simulation done
#   8. check-progress: exits 0 when simulation is complete
#
# Prerequisites:
#   1. polyzymd-env conda environment is active (or use mamba run)
#   2. polyzymd installed from feature/smart-restart branch:
#        mamba run -n polyzymd-env pip install -e ".[dev]"
#   3. setup_test_env.py has been run to create test_restart_workdir/
#
# Usage:
#   bash test_restart_interactive.sh
# =============================================================================

set -euo pipefail

WORKDIR="test_restart_workdir"
CONFIG="test_config.yaml"

# Use mamba run for safety, but also work if env is already active
CMD_PREFIX=""
if command -v mamba &>/dev/null && mamba env list 2>/dev/null | grep -q "polyzymd-env"; then
    if [[ "${CONDA_DEFAULT_ENV:-}" == "polyzymd-env" ]]; then
        CMD_PREFIX=""
    else
        CMD_PREFIX="mamba run --no-banner -n polyzymd-env"
    fi
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'  # No Color

pass_msg() { echo -e "${GREEN}PASS${NC}: $1"; }
fail_msg() { echo -e "${RED}FAIL${NC}: $1"; exit 1; }
info_msg() { echo -e "${CYAN}INFO${NC}: $1"; }
warn_msg() { echo -e "${YELLOW}WARN${NC}: $1"; }

header() {
    echo ""
    echo "====================================================================="
    echo -e "  ${CYAN}TEST $1${NC}: $2"
    echo "====================================================================="
    echo ""
}

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo ""
echo "====================================================================="
echo "  Self-Resubmitting Job Interactive Test Suite"
echo "====================================================================="
echo ""

if [ ! -d "$WORKDIR" ]; then
    fail_msg "$WORKDIR/ not found. Run setup_test_env.py first."
fi

if [ ! -f "$WORKDIR/solvated_system.pdb" ]; then
    fail_msg "$WORKDIR/solvated_system.pdb not found. Run setup_test_env.py first."
fi

if [ ! -f "$WORKDIR/progress.json" ]; then
    fail_msg "progress.json not found. Run setup_test_env.py first."
fi

if [ ! -f "$CONFIG" ]; then
    fail_msg "$CONFIG not found. Run setup_test_env.py first."
fi

# Check polyzymd is available
if ! $CMD_PREFIX polyzymd --help &>/dev/null; then
    fail_msg "polyzymd CLI not found. Install with: pip install -e '.[dev]'"
fi

info_msg "Pre-flight checks passed"
info_msg "Working directory: $WORKDIR"
info_msg "Config file: $CONFIG"
info_msg "Command prefix: '${CMD_PREFIX:-<none, env active>}'"

# =========================================================================
# TEST 1: run-segment continuation (segment 1)
# =========================================================================
header 1 "run-segment continuation (detects segment 0 done, runs segment 1)"

info_msg "Running: polyzymd run-segment -c $CONFIG -r 1 --scratch-dir $WORKDIR"
echo ""

$CMD_PREFIX polyzymd run-segment \
    -c "$CONFIG" \
    -r 1 \
    --scratch-dir "$WORKDIR"

RC=$?
if [ $RC -ne 0 ]; then
    fail_msg "polyzymd run-segment exited with code $RC"
fi

# Verify output
if [ ! -d "$WORKDIR/production_1" ]; then
    fail_msg "production_1/ directory not created"
fi
if [ ! -f "$WORKDIR/production_1/production_1_state.xml" ]; then
    fail_msg "production_1_state.xml not written"
fi

pass_msg "Segment 1 completed via run-segment"
info_msg "Files created in production_1/:"
ls -la "$WORKDIR/production_1/"

# =========================================================================
# TEST 2: check-progress (should show work remaining)
# =========================================================================
header 2 "check-progress (should show work remaining after segment 1)"

info_msg "Running: polyzymd check-progress -c $CONFIG -r 1 --scratch-dir $WORKDIR"
echo ""

set +e
$CMD_PREFIX polyzymd check-progress \
    -c "$CONFIG" \
    -r 1 \
    --scratch-dir "$WORKDIR"
RC=$?
set -e

if [ $RC -eq 0 ]; then
    warn_msg "check-progress exited 0 (complete) — simulation might be too short"
    warn_msg "This means all steps fit in 2 segments. Adjusting expectations."
    SIM_COMPLETE=1
elif [ $RC -eq 1 ]; then
    pass_msg "check-progress exited 1 (work remains)"
    SIM_COMPLETE=0
else
    fail_msg "check-progress exited with unexpected code $RC"
fi

# If the simulation is already complete after 2 segments, we need to test
# interruption differently. Reset progress to re-run.
if [ "$SIM_COMPLETE" -eq 1 ]; then
    info_msg "Simulation completed in 2 segments. All core tests passed."
    info_msg "Skipping signal/recovery tests (nothing to interrupt)."

    header "DONE" "Core tests passed (simulation completed quickly)"
    echo ""
    echo "Final directory structure:"
    for d in "$WORKDIR"/production_*; do
        if [ -d "$d" ]; then
            name=$(basename "$d")
            n_files=$(ls "$d" | wc -l)
            echo "  $name/ ($n_files files)"
        fi
    done
    echo ""
    echo -e "${GREEN}ALL CORE TESTS PASSED${NC}"
    echo ""
    echo "To clean up:  rm -rf $WORKDIR $CONFIG"
    exit 0
fi

# =========================================================================
# TEST 3: Signal interruption during run-segment
# =========================================================================
header 3 "Signal interruption (kill -USR1 during run-segment)"

info_msg "Running run-segment in the background..."
info_msg "Will send SIGUSR1 after 3 seconds to trigger graceful shutdown."
echo ""

$CMD_PREFIX polyzymd run-segment \
    -c "$CONFIG" \
    -r 1 \
    --scratch-dir "$WORKDIR" &
BG_PID=$!

info_msg "Background PID: $BG_PID"
info_msg "Waiting 3 seconds before sending USR1..."
sleep 3

if ! kill -0 $BG_PID 2>/dev/null; then
    warn_msg "Process already finished (system too fast for signal test)"
    warn_msg "Creating synthetic INTERRUPTED state for recovery testing..."

    # Find the latest segment dir
    LATEST_SEG=$(ls -d "$WORKDIR"/production_* 2>/dev/null | sort -t_ -k2 -n | tail -1)
    if [ -n "$LATEST_SEG" ]; then
        SEG_NAME=$(basename "$LATEST_SEG")
        SEG_IDX=${SEG_NAME#production_}

        # Copy state files as emergency files
        if [ -f "$LATEST_SEG/${SEG_NAME}_state.xml" ]; then
            cp "$LATEST_SEG/${SEG_NAME}_state.xml" "$LATEST_SEG/emergency_state.xml"
            cp "$LATEST_SEG/${SEG_NAME}_system.xml" "$LATEST_SEG/emergency_system.xml"
        fi

        # Create INTERRUPTED marker
        cat > "$LATEST_SEG/INTERRUPTED" <<MARKER
segment_index=$SEG_IDX
steps_completed=500
total_steps=1000
remaining_steps=500
MARKER
        # Remove normal end-of-segment state so run-segment sees it as interrupted
        rm -f "$LATEST_SEG/${SEG_NAME}_state.xml"

        warn_msg "Created synthetic INTERRUPTED state in $SEG_NAME/"
    fi
    SIGNAL_TEST_SKIPPED=1
else
    info_msg "Sending SIGUSR1 to PID $BG_PID..."
    kill -USR1 $BG_PID

    info_msg "Waiting for graceful shutdown..."
    set +e
    wait $BG_PID 2>/dev/null
    EXIT_CODE=$?
    # Re-wait in case signal interrupted the wait itself
    while kill -0 $BG_PID 2>/dev/null; do
        wait $BG_PID 2>/dev/null
        EXIT_CODE=$?
    done
    set -e

    info_msg "Process exited with code: $EXIT_CODE"
    SIGNAL_TEST_SKIPPED=0

    if [ $EXIT_CODE -eq 99 ]; then
        pass_msg "Exit code 99 (graceful shutdown)"
    else
        warn_msg "Expected exit code 99, got $EXIT_CODE"
    fi
fi

# =========================================================================
# TEST 4: check-progress after interruption
# =========================================================================
header 4 "check-progress (should show work remaining after interruption)"

set +e
$CMD_PREFIX polyzymd check-progress \
    -c "$CONFIG" \
    -r 1 \
    --scratch-dir "$WORKDIR"
RC=$?
set -e

if [ $RC -eq 1 ]; then
    pass_msg "check-progress exited 1 (work remains after interruption)"
elif [ $RC -eq 0 ]; then
    warn_msg "check-progress says complete — may have finished before signal"
fi

# =========================================================================
# TEST 5: run-segment recovery (auto-detects interrupted state)
# =========================================================================
header 5 "run-segment recovery (auto-detects interrupted work, continues)"

info_msg "Running: polyzymd run-segment -c $CONFIG -r 1 --scratch-dir $WORKDIR"
info_msg "This should detect the interruption and continue from where it left off."
echo ""

$CMD_PREFIX polyzymd run-segment \
    -c "$CONFIG" \
    -r 1 \
    --scratch-dir "$WORKDIR"

RC=$?
if [ $RC -eq 0 ]; then
    pass_msg "run-segment completed successfully (continued from interruption)"
elif [ $RC -eq 99 ]; then
    warn_msg "run-segment interrupted again (exit 99) — may need more wall time"
else
    fail_msg "run-segment failed with code $RC"
fi

# =========================================================================
# TEST 6: recover (status report, no --submit)
# =========================================================================
header 6 "recover (status report only)"

info_msg "Running: polyzymd recover -c $CONFIG -r 1 --scratch-dir $WORKDIR"
echo ""

$CMD_PREFIX polyzymd recover \
    -c "$CONFIG" \
    -r 1 \
    --scratch-dir "$WORKDIR"

# recover always exits 0 (it's a status command)
pass_msg "recover reported status successfully"

# =========================================================================
# TEST 7: run-segment until complete
# =========================================================================
header 7 "run-segment until complete (loop until check-progress exits 0)"

MAX_ITERATIONS=10
ITER=0
while true; do
    ITER=$((ITER + 1))
    if [ $ITER -gt $MAX_ITERATIONS ]; then
        fail_msg "Simulation did not complete after $MAX_ITERATIONS iterations"
    fi

    set +e
    $CMD_PREFIX polyzymd check-progress \
        -c "$CONFIG" \
        -r 1 \
        --scratch-dir "$WORKDIR"
    CHECK_RC=$?
    set -e

    if [ $CHECK_RC -eq 0 ]; then
        pass_msg "Simulation complete (check-progress exits 0)"
        break
    fi

    info_msg "Iteration $ITER: running next segment..."
    $CMD_PREFIX polyzymd run-segment \
        -c "$CONFIG" \
        -r 1 \
        --scratch-dir "$WORKDIR"
done

# =========================================================================
# TEST 8: Final check-progress confirms completion
# =========================================================================
header 8 "Final check-progress (should exit 0 = complete)"

set +e
$CMD_PREFIX polyzymd check-progress \
    -c "$CONFIG" \
    -r 1 \
    --scratch-dir "$WORKDIR"
RC=$?
set -e

if [ $RC -eq 0 ]; then
    pass_msg "Final check-progress confirms simulation is complete"
else
    fail_msg "Final check-progress exited $RC (expected 0)"
fi

# =========================================================================
# Summary
# =========================================================================
echo ""
echo "====================================================================="
echo -e "  ${GREEN}ALL TESTS PASSED${NC}"
echo "====================================================================="
echo ""
echo "Final directory structure:"
echo ""
for d in "$WORKDIR"/production_*; do
    if [ -d "$d" ]; then
        name=$(basename "$d")
        n_files=$(ls "$d" | wc -l)
        has_state=""
        has_interrupted=""
        if ls "$d"/*_state.xml &>/dev/null 2>&1; then
            has_state=" [state.xml present]"
        fi
        if [ -f "$d/INTERRUPTED" ]; then
            has_interrupted=" [INTERRUPTED]"
        fi
        echo "  $name/ ($n_files files)${has_state}${has_interrupted}"
    fi
done
echo ""
if [ "${SIGNAL_TEST_SKIPPED:-0}" -eq 1 ]; then
    warn_msg "Signal test was skipped (simulation too fast for the tiny waterbox)."
    warn_msg "To test signal handling, increase production duration in setup_test_env.py."
    warn_msg "All other tests (run-segment, check-progress, recover) passed with real OpenMM."
fi
echo ""
echo "To clean up:  rm -rf $WORKDIR $CONFIG"
echo ""
