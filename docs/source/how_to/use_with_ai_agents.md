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

## Refresh the files after you update PolyzyMD

Each file holds a marker line with the PolyzyMD version that wrote it and the
SHA-256 of the rest of the file, for example:

```
<!-- polyzymd-managed: version=1.4.0 sha256=8c3cf6f5... polyzymd agent sync rewrites this file unless you edit it -->
```

After you update PolyzyMD, rewrite the files to the installed version from
anywhere inside the project:

```bash
polyzymd agent sync --dry-run    # say what would change
polyzymd agent sync
```

```
unchanged: AGENTS.md
unchanged: CLAUDE.md
updated: .claude/skills/polyzymd-analyze/SKILL.md
updated: .agents/skills/polyzymd-analyze/SKILL.md
created: .claude/skills/polyzymd-organize/SKILL.md
created: .agents/skills/polyzymd-organize/SKILL.md
next: review and commit them in /home/me/Paper_1: git add ...
```

Sync commits nothing. Review the changes and commit them.

## Keep your own edits

You can add your own rules to `AGENTS.md`. Sync sees that the file no longer
matches the hash in its marker and leaves it alone:

```
skipped: AGENTS.md (edited; --force to overwrite)
```

`--force` overwrites edited files with the installed version, so your edits
are lost unless git has them. A file without a marker, such as an `AGENTS.md`
you wrote yourself, is skipped in the same way.

## Add the files to an existing project

A project or study made with an older PolyzyMD has none of these files. Run
`polyzymd agent sync` inside it to create them, then commit them.

See also: {doc}`project`, {doc}`study_folder` and the
{doc}`../reference/cli_reference` entry for `polyzymd agent sync`.
