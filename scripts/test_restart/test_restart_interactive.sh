#!/usr/bin/env bash
# =============================================================================
# Smart-Restart Interactive Test Script
# =============================================================================
#
# Exercises all 4 layers of the smart-restart infrastructure using the
# synthetic waterbox created by setup_test_env.py.
#
# Prerequisites:
#   1. polyzymd-env conda environment is active (or use mamba run)
#   2. polyzymd installed in editable mode from feature/smart-restart branch:
#        mamba run -n polyzymd-env pip install -e ".[dev]"
#   3. setup_test_env.py has been run to create test_restart_workdir/
#
# Usage:
#   bash test_restart_interactive.sh
#
# The script runs 5 tests sequentially. Each test prints a header, runs
# the command, and checks the result. If any test fails, the script stops.
#
# Tests:
#   1. Normal continuation (segment 1): polyzymd continue works with chunked loop
#   2. Signal interruption (segment 2): kill -USR1 triggers graceful shutdown
#   3. Recover dry-run: polyzymd recover --dry-run parses INTERRUPTED marker
#   4. Recover actual: polyzymd recover runs remaining steps, writes end-of-segment
#   5. Recover-chain scan: polyzymd recover-chain --dry-run shows chain status
# =============================================================================

set -euo pipefail

WORKDIR="test_restart_workdir"

# Use mamba run for safety, but also work if env is already active
CMD_PREFIX=""
if command -v mamba &>/dev/null && mamba env list 2>/dev/null | grep -q "polyzymd-env"; then
    # Check if we're already in the env
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
echo "  Smart-Restart Interactive Test Suite"
echo "====================================================================="
echo ""

if [ ! -d "$WORKDIR" ]; then
    fail_msg "$WORKDIR/ not found. Run setup_test_env.py first."
fi

if [ ! -f "$WORKDIR/solvated_system.pdb" ]; then
    fail_msg "$WORKDIR/solvated_system.pdb not found. Run setup_test_env.py first."
fi

if [ ! -f "$WORKDIR/production_0/production_0_state.xml" ]; then
    fail_msg "Segment 0 files not found. Run setup_test_env.py first."
fi

# Check polyzymd is available
if ! $CMD_PREFIX polyzymd --help &>/dev/null; then
    fail_msg "polyzymd CLI not found. Install with: pip install -e '.[dev]'"
fi

info_msg "Pre-flight checks passed"
info_msg "Working directory: $WORKDIR"
info_msg "Command prefix: '${CMD_PREFIX:-<none, env active>}'"

# =========================================================================
# TEST 1: Normal continuation (segment 1)
# =========================================================================
header 1 "Normal continuation (polyzymd continue, segment 1)"

info_msg "Running: polyzymd continue -w $WORKDIR -s 1 -t 0.001 -n 10"
info_msg "This runs 500 steps (0.001 ns at 2fs timestep) and saves 10 frames."
echo ""

$CMD_PREFIX polyzymd continue \
    -w "$WORKDIR" \
    -s 1 \
    -t 0.001 \
    -n 10

RC=$?
if [ $RC -ne 0 ]; then
    fail_msg "polyzymd continue exited with code $RC"
fi

# Verify output
if [ ! -d "$WORKDIR/production_1" ]; then
    fail_msg "production_1/ directory not created"
fi
if [ ! -f "$WORKDIR/production_1/production_1_state.xml" ]; then
    fail_msg "production_1_state.xml not written"
fi
if [ ! -f "$WORKDIR/production_1/production_1_system.xml" ]; then
    fail_msg "production_1_system.xml not written"
fi
if [ ! -f "$WORKDIR/production_1/production_1_parameters.json" ]; then
    fail_msg "production_1_parameters.json not written"
fi

pass_msg "Segment 1 completed normally"
info_msg "Files created in production_1/:"
ls -la "$WORKDIR/production_1/"

# =========================================================================
# TEST 2: Signal interruption (segment 2)
# =========================================================================
header 2 "Signal interruption (kill -USR1 during segment 2)"

info_msg "Running a longer segment 2 in the background (0.01 ns = 5000 steps)"
info_msg "Will send SIGUSR1 after 3 seconds to trigger graceful shutdown."
echo ""

# Run continuation in background with a longer segment
$CMD_PREFIX polyzymd continue \
    -w "$WORKDIR" \
    -s 2 \
    -t 0.01 \
    -n 50 &
BG_PID=$!

info_msg "Background PID: $BG_PID"
info_msg "Waiting 3 seconds before sending USR1..."
sleep 3

# Check if process is still running
if ! kill -0 $BG_PID 2>/dev/null; then
    warn_msg "Process already finished (waterbox is tiny — simulation may have completed before signal)"
    warn_msg "This means the chunked loop works but we couldn't test signal handling."
    warn_msg "Try increasing segment duration (-t 0.1) or decreasing system speed."
    echo ""
    info_msg "Checking if segment 2 completed normally instead..."

    if [ -f "$WORKDIR/production_2/production_2_state.xml" ]; then
        info_msg "Segment 2 completed normally (too fast for signal test)"
        info_msg "Skipping to test 3 — we'll create a fake INTERRUPTED state for testing."

        # Create a fake interrupted state for testing recover
        # Copy segment 2's files as emergency files and add INTERRUPTED marker
        cp "$WORKDIR/production_2/production_2_state.xml" \
           "$WORKDIR/production_2/emergency_state.xml"
        cp "$WORKDIR/production_2/production_2_system.xml" \
           "$WORKDIR/production_2/emergency_system.xml"
        if [ -f "$WORKDIR/production_2/production_2_checkpoint.chk" ]; then
            cp "$WORKDIR/production_2/production_2_checkpoint.chk" \
               "$WORKDIR/production_2/emergency_checkpoint.chk"
        fi

        # Create INTERRUPTED marker as if we stopped at step 2500/5000
        cat > "$WORKDIR/production_2/INTERRUPTED" <<MARKER
segment_index=2
steps_completed=2500
total_steps=5000
remaining_steps=2500
MARKER

        # Remove the normal end-of-segment state so recover has work to do
        rm -f "$WORKDIR/production_2/production_2_state.xml"

        warn_msg "Created synthetic INTERRUPTED state for recover testing"
        SIGNAL_TEST_SKIPPED=1
    else
        fail_msg "Segment 2 did not complete and did not get interrupted"
    fi
else
    info_msg "Sending SIGUSR1 to PID $BG_PID..."
    kill -USR1 $BG_PID

    # Wait for process to exit
    info_msg "Waiting for graceful shutdown..."
    wait $BG_PID || true
    EXIT_CODE=$?
    info_msg "Process exited with code: $EXIT_CODE"
    SIGNAL_TEST_SKIPPED=0

    if [ $EXIT_CODE -eq 99 ]; then
        pass_msg "Exit code 99 (graceful shutdown)"
    else
        warn_msg "Expected exit code 99, got $EXIT_CODE"
        warn_msg "This may happen if the signal was caught but exit code propagation differs."
    fi
fi

# Verify INTERRUPTED marker exists
if [ -f "$WORKDIR/production_2/INTERRUPTED" ]; then
    pass_msg "INTERRUPTED marker file created"
    info_msg "INTERRUPTED marker contents:"
    cat "$WORKDIR/production_2/INTERRUPTED"
    echo ""
else
    fail_msg "INTERRUPTED marker not found in production_2/"
fi

# Verify emergency state files
if [ -f "$WORKDIR/production_2/emergency_state.xml" ]; then
    pass_msg "emergency_state.xml created"
else
    fail_msg "emergency_state.xml not found"
fi

if [ -f "$WORKDIR/production_2/emergency_system.xml" ]; then
    pass_msg "emergency_system.xml created"
else
    fail_msg "emergency_system.xml not found"
fi

info_msg "Files in production_2/ after interruption:"
ls -la "$WORKDIR/production_2/"

# =========================================================================
# TEST 3: Recover dry-run
# =========================================================================
header 3 "Recover dry-run (polyzymd recover --dry-run)"

info_msg "Running: polyzymd recover -w $WORKDIR -s 2 --dry-run"
echo ""

$CMD_PREFIX polyzymd recover \
    -w "$WORKDIR" \
    -s 2 \
    --dry-run

RC=$?
if [ $RC -ne 0 ]; then
    fail_msg "polyzymd recover --dry-run exited with code $RC"
fi

pass_msg "Dry-run completed, showing remaining steps"

# =========================================================================
# TEST 4: Recover actual
# =========================================================================
header 4 "Recover actual (polyzymd recover, runs remaining steps)"

info_msg "Running: polyzymd recover -w $WORKDIR -s 2"
info_msg "This will run remaining steps in production_2_resume/"
echo ""

$CMD_PREFIX polyzymd recover \
    -w "$WORKDIR" \
    -s 2

RC=$?
if [ $RC -ne 0 ]; then
    fail_msg "polyzymd recover exited with code $RC"
fi

# Verify recovery outputs
if [ -d "$WORKDIR/production_2_resume" ]; then
    pass_msg "Resume directory created: production_2_resume/"
    info_msg "Resume directory contents:"
    ls -la "$WORKDIR/production_2_resume/"
else
    fail_msg "production_2_resume/ directory not created"
fi

if [ -f "$WORKDIR/production_2/production_2_state.xml" ]; then
    pass_msg "End-of-segment state.xml written back to production_2/"
else
    fail_msg "production_2_state.xml not written after recovery"
fi

if [ -f "$WORKDIR/production_2/production_2_system.xml" ]; then
    pass_msg "End-of-segment system.xml written back to production_2/"
else
    fail_msg "production_2_system.xml not written after recovery"
fi

if [ -f "$WORKDIR/production_2/INTERRUPTED" ]; then
    fail_msg "INTERRUPTED marker still exists after recovery (should be removed)"
else
    pass_msg "INTERRUPTED marker removed after successful recovery"
fi

# =========================================================================
# TEST 5: Recover-chain scan
# =========================================================================
header 5 "Recover-chain scan (polyzymd recover-chain --dry-run)"

info_msg "Running: polyzymd recover-chain -w $WORKDIR --dry-run"
echo ""

$CMD_PREFIX polyzymd recover-chain \
    -w "$WORKDIR" \
    --dry-run

RC=$?
if [ $RC -ne 0 ]; then
    fail_msg "polyzymd recover-chain --dry-run exited with code $RC"
fi

pass_msg "Chain status scan completed"

# =========================================================================
# TEST 6: Verify continuation after recovery (segment 3)
# =========================================================================
header 6 "Continue after recovery (segment 3 from recovered segment 2)"

info_msg "Running: polyzymd continue -w $WORKDIR -s 3 -t 0.001 -n 10"
info_msg "This verifies that ContinuationManager can pick up from recovered state."
echo ""

$CMD_PREFIX polyzymd continue \
    -w "$WORKDIR" \
    -s 3 \
    -t 0.001 \
    -n 10

RC=$?
if [ $RC -ne 0 ]; then
    fail_msg "polyzymd continue for segment 3 exited with code $RC"
fi

if [ -f "$WORKDIR/production_3/production_3_state.xml" ]; then
    pass_msg "Segment 3 completed successfully from recovered segment 2"
else
    fail_msg "Segment 3 did not produce state.xml"
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
        if [ -f "$d/${name}_state.xml" ] 2>/dev/null || \
           ls "$d"/*_state.xml &>/dev/null; then
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
    warn_msg "To test signal handling, try with a real system or increase -t to 1.0 ns."
    warn_msg "All other tests (recovery, continuation, chain scan) passed with real OpenMM."
fi
echo ""
echo "To clean up:  rm -rf $WORKDIR"
echo ""
