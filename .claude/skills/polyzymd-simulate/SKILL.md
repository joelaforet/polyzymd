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
user's PDB in `<study>/conditions/<name>/structures/` and set `pdb_path` on
the `type: protein` entry of `solutes:`.

**The scientific choices are the user's.** Every solute names its
`force_field` files and its `charges` (with the NAGL `model` file for
`method: nagl`), and `solvent` names its `force_field` and `water.model`.
PolyzyMD fills in none of them. Take them from the user or from a config of
theirs that names them, and say in your reply which ones the config uses. A
1.3.0rc4 config that names no force field or NAGL model was built with
`openff-2.0.0.offxml` and `openff-gnn-am1bcc-0.1.0-rc.3.pt`: keep those to
reproduce it, or ask the user before moving to newer ones. Never change one
to make a build pass.

A config with `enzyme:`, `substrate:`, `polymers:` or `force_field:` is from
1.3.0rc4 or earlier and is refused: rewrite it with
`docs/source/how_to/rewrite_config_rc5.md`.

Then, with `C=trpcage/conditions/water/config.yaml`:

```bash
polyzymd validate -c $C          # check the config; names every wrong key and missing file
polyzymd build -c $C --dry-run   # what will be built: components, counts, seeds, folders
polyzymd build -c $C -r 1        # build replicate 1 and run the physics audit
polyzymd run -c $C -r 1          # run replicate 1 here (engine from the config)
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

- **Protein:** `pdb_path` of the `type: protein` solute, with
  `charges: {method: library}`. Prepare a raw PDB entry in this order:
  1. Keep the protein you mean to simulate. If the file holds several copies
     of it, keep one copy (or the assembly you want).
  2. List the residue names of the `HETATM` records:
     `grep "^HETATM" raw.pdb | cut -c18-20 | sort | uniq -c`. Delete the
     waters (`HOH`), ions, crystallization additives and ligands. A ligand
     you want goes in a `type: small_molecule` solute as its own SDF. If nothing else is
     listed, `grep -v "^HETATM" raw.pdb > protein.pdb` does this step.
     Modified amino acids are also `HETATM` records. PolyzyMD cannot yet
     parameterize them (issue #228), so keep them and decide each one on
     purpose:
     - Selenomethionine (`MSE`) is usually a phasing substitution for the
       native methionine. Convert it to `MET` by hand before clean-pdb,
       which would drop the selenium without adding a sulfur: change the
       records to `ATOM`, the residue name `MSE` to `MET`, the atom name
       `SE` to `SD` and its element `SE` to `S`.
     - clean-pdb replaces any other modified residue with its parent
       residue (phosphoserine `SEP` becomes serine), which changes the
       chemistry. Tell the user before you simulate it that way.
  3. Run `polyzymd clean-pdb -i protein.pdb -o clean.pdb`. It replaces
     nonstandard residues and adds hydrogens (pH 7.4 unless you give `--ph`).
     It does not model missing residues or atoms; if `REMARK 465` or
     `REMARK 470` lists any, see
     `docs/source/how_to/troubleshoot_openff_pdb_ingestion.md`. For a chain
     that only ends early, follow "Add the terminal oxygen of a chain
     that ends early" in
     `docs/source/tutorials/prepare_pdb_for_openff.md`.
  4. Keep the chain IDs as they are, and a `TER` record between chains and
     between disconnected fragments. Do not merge chains into one or give
     them one chain ID by editing the PDB: a merged chain gets peptide bonds
     that do not exist, and chains joined by disulfides need distinct IDs
     for the build to number them apart. Each chain, or each group of chains
     joined by disulfide bonds, loads as one protein molecule; the build
     keeps every protein molecule, puts them all on chain A and numbers
     their residues on from one molecule to the next (see
     `docs/source/explanation/residue_assignment.md`). Take residue numbers
     for restraints and analysis selections from the built
     `solvated_system.pdb`, not from the raw file.
- **Ligand:** a `type: small_molecule` entry of `solutes:` with a docked,
  protonated SDF (`sdf_path`, `residue_name`, `force_field`, and `charges`
  such as `{method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}`). One ligand
  only.
- **Co-solvent from a SMILES:** under `solvent.co_solvents`, give `name`,
  `smiles`, `residue_name` (3 letters) and one of `count`, `concentration`
  (M) or `mole_fraction`, and `charges`: `{method: preset}` for a co-solvent
  PolyzyMD ships (dmso, urea, ...), otherwise `nagl` with its `model`, or
  `am1bcc` (needs AmberTools). The build warns you to check NAGL charges. Write a charged molecule with its counter-ion
  (`CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]`), or without it and keep
  `solvent.ions.neutralize: true`: either way the system is neutral.
- **Polymers:** a `type: polymer` entry; see the polymers how-to. Cached
  chains take `charges: {method: preset}`.
- **Restraints:** each atom names the solute it belongs to (`component:`),
  and the selection is matched inside that solute only.
- **Engine:** `engine: gromacs` or `engine: openmm`. The page "GROMACS and
  OpenMM" in the reference says how each setting runs on each engine and
  which mappings are approximate.

## 3. Read the build record

Every build ends with the physics audit: the built system is compared, term
by term, with what the config names. Open `build_provenance.md` in the
replicate folder and tell the user what it says: the components with their
force fields and charges, and any WARN line.

If the audit fails, the build stops and names each failing check
(`FAIL <check> (<component>): <detail>`); the full record is in
`build_provenance.failed.md` (GROMACS: `build_provenance.md`, marked
failed). Most failures mean the config does not say
what the user meant: fix the config, not the check. A `modeling_exceptions`
entry accepts a failure on purpose, and its `reason` is published with the
config. Add one only when the user has told you the departure is deliberate,
and write their reason, not yours. Check IDs are listed in
`docs/source/reference/physics_checks.md`.

Never edit an exported `.top`, `.itp`, `.mdp` or `system.xml` by hand: runs
refuse files the audit did not vouch for. Change the config and build again.

## 4. Run on a cluster

```bash
polyzymd submit -c $C -r 1-5 --dry-run   # show what would be submitted first
polyzymd submit -c $C -r 1-5             # submit; jobs resubmit themselves
polyzymd status -c $C
```

Check `--preset` and its SLURM account before submitting.

## 5. Analyse

List the analyses under `analyses:` in `project.yaml` (such as `rg: {}`),
commit, and use the `polyzymd-analyze` skill:

```bash
git add -A && git commit -m "Conditions and analyses"
polyzymd analyze --project .
```

## 6. When something fails

- Read the `error:` and `fix:` lines; they name the key or file.
- `polyzymd validate` refuses unknown keys: fix the spelling, or remove the
  key.
- `run` or `submit` says no passing physics audit vouches for the files: the
  build is missing, older than 1.3.0rc5, stale for the current config, or a
  file was edited. Do what its `fix:` line says (rebuild, or restore the
  file).
- `polyzymd status` shows `SYSTEM DIVERGED FROM THE AUDITED BUILD`: the
  run integrated a system other than the audited one. Tell the user; do not
  resubmit or rebuild on your own.
- `grompp` stops on any warning: read it. A net charge means the system is
  not neutral.
- A topology MDAnalysis cannot read: `polyzymd analysis-topology --overwrite RUN_DIR`.
