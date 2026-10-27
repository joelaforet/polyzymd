#!/usr/bin/env bash
# PreToolUse gate on Bash. Blocks `git commit` until ruff, black and the
# scientific analysis tests pass. See .claude/settings.json for the wiring and
# .claude/skills/livecoms-check/SKILL.md for what the tests are protecting.
#
# Exit 0 lets the command through. Exit 2 blocks it and shows stderr to Claude.

set -uo pipefail

payload=$(cat)

block() {
    echo "commit gate: $1" >&2
    echo "Fix it and commit again. Run the livecoms-check skill if the failure is scientific." >&2
    exit 2
}

# Pull tool_input.command out of the hook payload. A payload we cannot parse is
# treated as a commit, so a malformed or unexpected payload blocks rather than
# waves the commit through. The no-Python path does the same by matching against
# the raw payload text.
command_text=""
if python_bin=$(command -v python3 || command -v python); then
    if ! command_text=$("$python_bin" -c '
import json, sys
payload = json.load(sys.stdin)
print(payload.get("tool_input", {}).get("command", ""))
' <<<"$payload" 2>/dev/null); then
        block "could not parse the hook payload, so the commit is not allowed through unchecked."
    fi
else
    command_text="$payload"
fi

# Match git commit in the forms that actually reach the shell, including
# `git -C <path> commit` and `git --no-pager commit`, at the start of the
# command or after a shell separator.
git_commit_re='(^|[;&|]\s*)git(\s+-[A-Za-z-]+(\s+[^[:space:]]+)?)*\s+commit\b'
if ! grep -Eq "$git_commit_re" <<<"$command_text"; then
    exit 0
fi

repo_root="${CLAUDE_PROJECT_DIR:-}"
if [[ -z "$repo_root" ]]; then
    repo_root=$(git rev-parse --show-toplevel 2>/dev/null) || {
        echo "commit gate: not inside a git repository, cannot locate the pixi environments." >&2
        exit 2
    }
fi

# Worktrees have no .pixi of their own. They share the main checkout's
# environments, which sit next to the common git directory.
pixi_root="$repo_root"
if [[ ! -d "$pixi_root/.pixi/envs" ]]; then
    common_dir=$(git -C "$repo_root" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)
    if [[ -n "$common_dir" ]]; then
        pixi_root=$(dirname "$common_dir")
    fi
fi

build_bin="$pixi_root/.pixi/envs/build/bin"
test_python="$pixi_root/.pixi/envs/test/bin/python"

for tool in "$build_bin/ruff" "$build_bin/black"; do
    [[ -x "$tool" ]] || block "cannot find $tool. Check that the build environment is installed."
done

cd "$repo_root" || block "cannot enter $repo_root."

if ! output=$("$build_bin/ruff" check src tests 2>&1); then
    block "ruff check failed."$'\n'"$output"
fi

if ! output=$("$build_bin/black" --check src tests 2>&1); then
    block "black --check failed."$'\n'"$output"
fi

if [[ -d tests/analyses/scientific ]]; then
    [[ -x "$test_python" ]] || block "cannot find $test_python. Check that the test environment is installed."
    if ! output=$(PYTHONPATH="$repo_root/src" "$test_python" -m pytest tests/analyses/scientific -q 2>&1); then
        block "the scientific analysis tests failed."$'\n'"$output"
    fi
fi

exit 0
