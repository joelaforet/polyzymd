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

## The 1.3 release work

`analyses_refactor` is the trunk of the 1.3 work. Pull requests on it form a
stack: `analyses/project-studies` (#163) targets `analyses_refactor`, and each
later branch (`analyses/contacts-zero-control` and the branches stacked on
it) targets the branch before it. `feature/v1.3.0-rc5` and `main` have commits that
`analyses_refactor` lacks; Joe syncs the trunks before the 1.3.0 tag.

- Branch from the top of the stack, or from the branch the maintainer names,
  as `analyses/<item>`. Open the pull request as a draft against that branch,
  never against `feature/v1.3.0-rc5` or `main`.
- Never commit on `analyses_refactor` itself. Joe merges.
- Never put plan, wave or slice names or finding IDs in files. Name tests
  and CHANGELOG entries by the behavior they check, and put each test in the
  module that mirrors the source module it checks (`tests/` mirrors
  `src/polyzymd/`).
