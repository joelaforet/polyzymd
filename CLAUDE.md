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

Do not run `pixi install`. The environments are already solved and a reinstall
costs an hour.

## Worktrees

A git worktree created from this repository shares the main checkout's pixi
environments, and those environments contain an editable install that points at
`~/Shirts-Lab-Linux/polyzymd/src`. Importing `polyzymd` from a worktree without
help gives you the main checkout's code, not yours. Set `PYTHONPATH` to your own
`src` and call the environment's interpreter directly.

```bash
# from inside the worktree
PYTHONPATH=$PWD/src ~/Shirts-Lab-Linux/polyzymd/.pixi/envs/test/bin/python -m pytest tests/analyses -q
PYTHONPATH=$PWD/src ~/Shirts-Lab-Linux/polyzymd/.pixi/envs/analysis/bin/python -c "..."
~/Shirts-Lab-Linux/polyzymd/.pixi/envs/build/bin/ruff check src tests
~/Shirts-Lab-Linux/polyzymd/.pixi/envs/build/bin/black --check src tests
```

Before trusting a test run from a worktree, print `polyzymd.__file__` and check
that it points inside the worktree.

## Skills

Repository skills live in `.claude/skills/`. Invoke `livecoms-check` before you
commit any change to `src/polyzymd/analyses/`. It holds the statistical rules
this project is held to, the fields every result must carry, and the known
answers the scientific tests check against.

## The analyses refactor

Work on `src/polyzymd/analyses/` follows the checklist in
`docs/planning/analyses_refactor.md`, which comes from the audit in
`docs/planning/analyses_audit_2026-09-11.md`. One checklist item is one branch,
one session and one pull request.

- Branch from `analyses_refactor` and name the branch `analyses/<item>`, using
  the name the checklist gives.
- Open the pull request against `analyses_refactor`, never against
  `feature/v1.3.0-rc5` or `main`.
- Never commit on `analyses_refactor` itself.
