@AGENTS.md

# Claude Code specifics

`AGENTS.md`, imported above, holds the repository rules that every coding agent
follows, including pixi environments and worktrees. This file adds the parts
that only apply to Claude Code.

## Skills

Repository skills live in `.claude/skills/`:

- `polyzymd-simulate`: write a config, validate, build, run, submit and check
  status.
- `polyzymd-analyze`: run analyses on a study and read the reports.
- `livecoms-check`: the statistical rules, the fields every result must carry
  and the known answers the scientific tests check. Invoke it before you
  commit any change to `src/polyzymd/analyses/`.

## The 1.4 release work

`release/1.4` is the trunk of the 1.4 work (covalently modified proteins on
the mBuild biopolymers fork). Branch from it, or from the branch the
maintainer names, and open the pull request as a draft against that branch,
never against `main`.

- Never commit on `release/1.4` itself. Joe merges.
- Never put plan, wave or slice names or finding IDs in files. Name tests
  and CHANGELOG entries by the behavior they check, and put each test in the
  module that mirrors the source module it checks (`tests/` mirrors
  `src/polyzymd/`).
