# Contributing to PolyzyMD

This guide is for developers and scientific contributors who want to make a
small, reviewable change to PolyzyMD and submit it cleanly.

## Start with the contributor environment

Set up the project as described in {doc}`setup`.

For most changes, your loop is:

1. create a branch
2. make one focused change
3. run the relevant checks
4. open a pull request with a clear explanation of why the change matters

## What good contributions look like

- one clear purpose per branch or PR
- changes that follow existing patterns instead of introducing a parallel style
- terminal colors from the shared helpers in `polyzymd/cli/colors.py`, never ad hoc ANSI escape codes
- docs and tests updated when behavior changes
- commands verified in a `pixi` environment

## Recommended verification commands

From the repository root:

```bash
pixi run -e build pytest tests/ -v
pixi run -e build ruff check src/
pixi run -e build make -C docs clean html
```

Run the narrowest useful subset for small changes, then run the broader set
before opening a PR if your change touches core workflows.

## A concrete example: add a new co-solvent

One common contribution is extending the built-in co-solvent library.

### 1. Add the library entry

Edit `src/polyzymd/data/cosolvent_library.py` and add a new item to
`COSOLVENT_LIBRARY`:

```python
"formamide": CoSolventData(
    name="Formamide",
    smiles="C(=O)N",
    density=1.133,
    molar_mass=45.04,
    residue_name="FMD",
    common_names=("methanamide",),
),
```

Use:

- a lowercase config name
- a valid SMILES string
- a literature or database-backed density
- a molar mass in g/mol
- a 3-letter residue name that is not an amino-acid, nucleic-acid, water or
  ion name (the PDB chemical-component code, where there is one)
- common alternative names that should resolve to the same solvent

### 2. Generate the solvent SDF

Run the generator from the repository root:

```bash
pixi run -e build python scripts/generate_solvent_sdfs.py
```

This creates a charged SDF in `src/polyzymd/data/solvents/`.

### 3. Test the new solvent in a config

Add it to a simple config:

```yaml
solvent:
  force_field: ["openff-2.3.0.offxml"]
  water:
    model: "tip3p"
  co_solvents:
    - name: "formamide"
      mole_fraction: 0.05
      charges: {method: preset}   # the charges in the generated SDF
```

Then verify the build path:

```bash
pixi run -e build polyzymd build -c test_config.yaml --dry-run
```

### 4. Include the right context in your PR

Call out:

- what you added
- the source of the density or other chemical constants
- any scientific caveats reviewers should know
- the commands you ran to verify the change

## Other common contribution areas

- docs cleanup and reorganization
- new analysis workflows or plotting support
- new comparison metrics
- HPC preset updates
- packaging and release automation

For extension-specific guidance, use the pages in
{doc}`index`.

## Writing documentation

The docs follow [Diataxis](https://diataxis.fr/): put a new page in the
section that matches what its reader is doing.

| Section | The reader wants to | A page here |
|---|---|---|
| Tutorials | learn, by doing a whole lesson with a known result | starts from nothing, says what it teaches and how long it takes |
| How-to guides | do one task they already understand | starts from the task, gives the steps, links to the reference |
| Reference | look up a fact | lists every key, command, file or check, and the evidence behind them |
| Explanation | understand why | discusses design and reasoning; has no steps to follow |

- Titles are in sentence case ("Run a simulation longer"); proper names keep
  their capitals (PolyzyMD, OpenMM, GROMACS, OpenFF, SLURM).
- Admonitions, one purpose each: `note` for context, `tip` for a shortcut or
  the environment to use ("Environment setup"), `important` for a step that
  must not be skipped, `warning` for anything that loses data or time.
- Add each new page to the toctree of its section's `index.md`, under the
  group (`:caption:`) it belongs to; the caption labels the group in the
  sidebar.
- Colours, the logo and the landing page live in `docs/source/_static/custom.css`
  and `docs/source/index.rst`; use the `--pz-*` colour variables rather than
  new colours.

## Before opening a PR

Check that your change:

- follows existing naming and file layout patterns
- does not commit generated junk or local cache files
- updates docs when user-facing behavior changes
- keeps experimental analysis workflows clearly labeled as experimental

## Related pages

- contributor environment: {doc}`setup`
- architecture overview: {doc}`architecture`
- packaging notes: {doc}`packaging`

<!-- IMAGE OPPORTUNITY: Add a contributor workflow figure showing `branch ->
edit -> test -> docs build -> PR`, with a small inset for the co-solvent data
path from library entry to generated SDF. -->
