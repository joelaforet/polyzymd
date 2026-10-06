---
name: polyzymd-simulate
description: Set up, build and run a molecular dynamics simulation with PolyzyMD, on OpenMM or GROMACS, locally or on SLURM. Use when a user wants to simulate a protein in water, with or without a ligand, co-solvents (a SMILES is enough) or polymers. Start from a project (polyzymd project init) and add conditions with polyzymd study add-condition; do not write a config from nothing.
---

# polyzymd-simulate: from a structure to trajectories

## 1. Make a project, then add conditions

Work top down: a project (one paper) holds studies (one protein each), and a
study holds conditions (one `config.yaml` each).

```bash
polyzymd project init my_paper --study trpcage    # project.yaml and the study folder trpcage/
cd my_paper
polyzymd study add-condition Water --new --study trpcage                 # template config to fill in
polyzymd study add-condition Water --config ~/polyzymd/examples/quickstart/config.yaml --study trpcage  # or copy a config
polyzymd study add-condition "Urea 2 M" --from Water --study trpcage     # copy a condition, then edit what differs
polyzymd project add-study calb --project .       # another protein
```

`examples/quickstart/config.yaml` is a protein in water with 0.15 M NaCl
(Trp-cage, OpenMM on the CPU; `config_gromacs.yaml` for GROMACS). Put the
user's PDB in `<study>/conditions/<name>/structures/` and set
`enzyme.pdb_path`. Then, with `C=trpcage/conditions/water/config.yaml`:

```bash
polyzymd validate -c $C          # check the config; names every wrong key and missing file
polyzymd build -c $C --dry-run   # what will be built: components, counts, seeds, folders
polyzymd run -c $C -r 1          # build and run replicate 1 here (engine from the config)
polyzymd status -c $C --no-slurm
```

Runs go into the git-ignored `runs/<study>/<name>/` of the project unless the
config sets `scratch_directory`. Tell the user plainly: trajectories can use
a lot of disk space, so on a cluster set `scratch_directory` to scratch
storage. `polyzymd init` is retired.

Paths in a config (`pdb_path`, `projects_directory`, ...) are relative to the
config's folder. `-r 1-3` runs replicates 1 to 3. The replicate number seeds
the starting structure (Packmol and polymer draws), the initial velocities
and the thermostat noise, so a replicate is reproducible on the same
hardware and software.

## 2. Change the system

- **Protein:** `enzyme.pdb_path`. Clean a raw PDB first with
  `polyzymd clean-pdb -i raw.pdb -o clean.pdb`. All protein residues must be
  in chain A.
- **Ligand:** a `substrate:` block with a docked, protonated SDF
  (`sdf_path`, `charge_method: nagl`, `residue_name`). One ligand only.
- **Co-solvent from a SMILES:** under `solvent.co_solvents`, give `name`,
  `smiles`, `residue_name` (3 letters) and one of `count`, `concentration`
  (M) or `mole_fraction`. Charges are NAGL unless you set
  `charge_method: am1bcc` (needs AmberTools); the build warns you to check
  them. Write a charged molecule with its counter-ion
  (`CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]`), or without it and keep
  `solvent.ions.neutralize: true`: either way the system is neutral.
- **Polymers:** see the polymers how-to; not needed for anything else.
- **Engine:** `engine: gromacs` or `engine: openmm`. The page "GROMACS and
  OpenMM" in the reference says how each setting runs on each engine and
  which mappings are approximate.

## 3. Run on a cluster

```bash
polyzymd submit -c $C -r 1-5 --dry-run   # show what would be submitted first
polyzymd submit -c $C -r 1-5             # submit; jobs resubmit themselves
polyzymd status -c $C
```

Check `--preset` and its SLURM account before submitting.

## 4. Analyse

List the analyses under `analyses:` in `project.yaml` (such as `rg: {}`),
commit, and use the `polyzymd-analyze` skill:

```bash
git add -A && git commit -m "Conditions and analyses"
polyzymd analyze --project .
```

## 5. When something fails

- Read the `error:` and `fix:` lines; they name the key or file.
- `polyzymd validate` refuses unknown keys: fix the spelling, or remove the
  key.
- `grompp` stops on any warning: read it. A net charge means the system is
  not neutral.
- A topology MDAnalysis cannot read: `polyzymd analysis-topology --overwrite RUN_DIR`.
