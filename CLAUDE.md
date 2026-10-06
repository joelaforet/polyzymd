@AGENTS.md

# Claude Code specifics

`AGENTS.md`, imported above, holds the repository rules that every coding agent
follows, including pixi environments and worktrees. This file adds the parts
that only apply to Claude Code.

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
