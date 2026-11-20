# Use PolyzyMD with AI coding agents

AI coding agents such as Claude Code and Codex read instructions from the
folder they work in. PolyzyMD writes these instructions into every new project
and study, so an agent that opens the folder knows the `polyzymd` commands and
the rules of the project before it writes any code.

## What PolyzyMD writes

`polyzymd project init` writes these files at the top of the project, and
`polyzymd study init` writes them at the top of a study that is not part of a
project:

| File | Read by | What it holds |
|------|---------|---------------|
| `AGENTS.md` | Codex, Claude Code (through `CLAUDE.md`) and other agents | The layout of the project, the core commands and the rules |
| `CLAUDE.md` | Claude Code | One line, `@AGENTS.md`, which loads `AGENTS.md` |
| `.claude/skills/<name>/SKILL.md` | Claude Code | The PolyzyMD skills, such as `polyzymd-analyze` and `polyzymd-simulate` |
| `.agents/skills/<name>/SKILL.md` | Codex and other agents that read `.agents/skills` | The same skills |

The files are part of the first commit. `AGENTS.md` tells the agent, among
other rules, to run `polyzymd analyze --list` before writing MDAnalysis code,
to put one-off scripts in `<study>/sandbox/`, never to loop over replicates
by hand, never to choose a force field, charge model, water model or
temperature for you, and to commit inputs after each decision with a message
that says why.

Freeze does not deposit these files, and changing them does not count as a
change to the inputs of an analysis.

## Start the agent in the project

Start the agent in the project folder, or in one of its studies, so it reads
`AGENTS.md`:

```bash
cd Paper_1
claude        # or: codex
```

If the project lives on a network mount (sshfs), the agent reads and edits
files through the mount but should run git and `polyzymd` on the machine that
hosts the files. `AGENTS.md` says so.

See also: {doc}`project` and {doc}`study_folder`.
