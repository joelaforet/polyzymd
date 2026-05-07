# Experimental analyses

PolyzyMD v1.3 ships the stable analysis plugins documented in the analysis
reference pages. The following experimental analyses were removed from the
shipped v1.3 code and are not available as active plugins:

- contacts binding preference
- exposure dynamics
- binding free energy
- polymer affinity
- polymer bridging

These archived features are retained only for historical review. They are not
maintained as part of the v1.3 analysis workflow, and this page intentionally
does not reproduce the old how-to content.

To inspect the archived implementation, fetch tags and switch to the archive
tag in a detached worktree state:

```bash
git fetch origin --tags
git switch --detach archive_experimental_analysis
```

If the historical branch is available locally or on your remote, you can also
inspect the last active development branch:

```bash
git switch feature/mda-analysis-migration
```

Return to your working branch before continuing v1.3 development.
