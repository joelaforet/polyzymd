---
name: polyzymd-simulate
description: Set up, build and run a molecular dynamics simulation with PolyzyMD, on OpenMM or GROMACS, locally or on SLURM. Use when a user wants to simulate a protein in water, with or without a ligand, co-solvents (a SMILES is enough) or polymers. Start from examples/quickstart/; do not write a config from nothing.
---

# polyzymd-simulate: from a structure to trajectories

## 1. Start from the quickstart example

`examples/quickstart/` is a protein in water with 0.15 M NaCl: Trp-cage,
`config.yaml` for OpenMM and `config_gromacs.yaml` for GROMACS. Copy the
folder, put your structure in place of `trpcage.pdb`, and edit the config.

```bash
polyzymd validate -c config.yaml          # check the config; names every wrong key
polyzymd build -c config.yaml --dry-run   # what will be built: components, counts, seeds
polyzymd run -c config.yaml -r 1          # build and run replicate 1 here (engine from the config)
polyzymd status -c config.yaml --no-slurm
```

Paths in a config (`pdb_path`, `projects_directory`, ...) are relative to the
config's folder. `-r 1-3` runs replicates 1 to 3. The replicate number seeds
the starting structure (Packmol and polymer draws); dynamics noise is drawn
afresh in each run.

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
polyzymd submit -c config.yaml -r 1-5 --dry-run   # show what would be submitted first
polyzymd submit -c config.yaml -r 1-5             # submit; jobs resubmit themselves
polyzymd status -c config.yaml
```

Check `--preset` and its SLURM account before submitting.

## 4. Analyse

Make a study of the conditions and use the `polyzymd-analyze` skill:

```bash
polyzymd study init my_study --condition "Water=water/config.yaml" --condition "SDS=sds/config.yaml"
polyzymd analyze rg --study my_study
```

## 5. When something fails

- Read the `error:` and `fix:` lines; they name the key or file.
- `polyzymd validate` refuses unknown keys: fix the spelling, or remove the
  key.
- `grompp` stops on any warning: read it. A net charge means the system is
  not neutral.
- A topology MDAnalysis cannot read: `polyzymd analysis-topology --overwrite RUN_DIR`.
