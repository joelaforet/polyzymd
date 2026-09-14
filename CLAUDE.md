@AGENTS.md

# Claude Code specifics

`AGENTS.md`, imported above, holds the repository rules that every coding agent
follows. This file adds the parts that only apply to Claude Code.

## Pixi environments

Three environments matter. The default environment has no numpy, so analysis
code and the test suite both fail there.

| Task | Environment | Example |
|------|-------------|---------|
| Run pytest | `test` | `pixi run -e test pytest tests/analyses -q` |
| Run analysis or library code | `analysis` | `pixi run -e analysis python -c "import polyzymd"` |
| Run ruff and black | `build` | `pixi run -e build ruff check src tests` |

Do not run `pixi install` or change `pixi.toml` during a session. Re-solving
the environments is slow and is the maintainer's job.

## Worktrees

A git worktree created from this repository shares the main checkout's pixi
environments, and those environments contain an editable install that points at
the main checkout's `src`. Importing `polyzymd` from a worktree without help
gives you the main checkout's code, not yours. Set `PYTHONPATH` to your own
`src` and call the environment's interpreter directly. The main checkout is the
parent of the common git directory:

```bash
# from inside the worktree
MAIN=$(dirname "$(git rev-parse --path-format=absolute --git-common-dir)")
PYTHONPATH=$PWD/src $MAIN/.pixi/envs/test/bin/python -m pytest tests/analyses -q
PYTHONPATH=$PWD/src $MAIN/.pixi/envs/analysis/bin/python -c "..."
$MAIN/.pixi/envs/build/bin/ruff check src tests
$MAIN/.pixi/envs/build/bin/black --check src tests
```

Before trusting a test run from a worktree, print `polyzymd.__file__` and check
that it points inside the worktree.

## Skills

Repository skills live in `.claude/skills/`. Invoke `livecoms-check` before you
commit any change to `src/polyzymd/analyses/`. It holds the statistical rules
this project is held to, the fields every result must carry, and the known
answers the scientific tests check against.

## The analyses refactor

Work on `src/polyzymd/analyses/` is split into items, each one branch, one
session and one pull request. The maintainer keeps the checklist outside the
repository and names the item when starting a session.

- Branch from `analyses_refactor` and name the branch `analyses/<item>`.
- Open the pull request against `analyses_refactor`, never against
  `feature/v1.3.0-rc5` or `main`.
- Never commit on `analyses_refactor` itself.
