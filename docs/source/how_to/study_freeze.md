# Publish a study with `polyzymd study freeze`

Freeze a {term}`study` when its analyses are done and you want to deposit it
with a paper. A reader can then make its figures from the deposit, and run
its analyses again from the trajectories. To create the folder, see
{doc}`study_folder`. For what the deposit contains and why, see
{doc}`../explanation/study_folders`.

A study that is part of a {term}`project` is frozen with the project. In that
case, use `polyzymd project freeze` (see {doc}`project`). `polyzymd study
freeze` refuses such a study.

:::{admonition} Environment Setup
:class: tip

The commands on this page assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

(study-metadata)=
## 1. Fill in the metadata

Add a `metadata:` block to `study.yaml`. These fields make the deposit
findable and citable:

```yaml
metadata:
  doi: "10.5281/zenodo.NNNNNNN"            # the DOI of the study, reserved in Zenodo (step 4)
  title: "LipA with SBMA-EGMA copolymers at 363 K"
  description: "Lipase A simulated without polymer and with SBMA-EGMA copolymers."
  purpose: "To test whether the copolymer stabilizes the lid at high temperature."
  keywords: [molecular dynamics, lipase, polymer, OpenMM]
  system_type: [protein, polymer]
  authors:
    - {family-names: Doe, given-names: Jane, orcid: "0000-0000-0000-0000",
       affiliation: "Example University"}
  license: {data: CC-BY-4.0, code: MIT}     # SPDX identifiers; match LICENSE-data and LICENSE-code
  funding:
    - {funder: "...", award: "...", funder_doi: "10.13039/..."}
  related:
    paper: {title: "...", status: in-preparation, doi: "10.XXXX/placeholder"}
    trajectories:                           # one entry for each trajectory deposit
      - {doi: "10.5281/zenodo.1234567", conditions: [No polymer, SBMA 50%]}
    experimental:
      - {doi: "10.1021/...", description: "Measured activity at 363 K"}
  zenodo: {communities: [], access_right: open}
```

| Field | Why |
|---|---|
| `doi` | The persistent identifier of the study (FAIR F1). Freeze writes it into `CITATION.cff`, `.zenodo.json` and the manifest. Reserve it in Zenodo before you publish |
| `title`, `description`, `keywords`, `authors` | Findability (FAIR F2) and citation |
| `purpose` | The main purpose of the simulations, part of the minimum metadata of Amaro et al. (2025) |
| `system_type` | The kinds of molecule simulated, as the Communications Biology checklist asks (3a) |
| `license` | How others can reuse the data and the code (FAIR R1.1) |
| `related.paper` | The DOI of the paper, when it has one. Until then, `status: in-preparation` |
| `related.trajectories` | The DOIs of the trajectory deposits, which you upload separately |
| `related.experimental` | The experiments that the simulations are compared with (checklist 2a) |

PolyzyMD refuses an unknown key, and gives the nearest spelling. A missing
field gives a warning and a `TODO` placeholder in the files. It does not stop
the freeze. So you can freeze early, and freeze again when you fill the gaps.
`polyzymd study check` prints the number of gaps that remain.

## 2. Prepare the study

1. Run every analysis:

   ```bash
   polyzymd analyze --study lipase_363K
   ```

2. If the replicates have no recorded trajectory hashes, record them once.
   Freeze warns when they are missing. See {doc}`study_folder`.

   ```bash
   polyzymd hash-trajectories --study lipase_363K
   ```

3. Put every file that you want to publish in a folder that freeze
   deposits. Freeze deposits only these names:

   - `study.yaml` (or `project.yaml`), `README*`, `LICENSE*` and the files
     that freeze writes;
   - the files under the folders that `study init` and `project init`
     make: `conditions/`, `structures/`, `analyses/`, `stats/`, `figures/`,
     `results/` and `environment/`. Reference structures go in
     `structures/`, and data files that a function reads go in
     `analyses/data/`.

   Freeze does not deposit other files, such as notes or a copied
   trajectory, and prints one `not deposited:` warning that names them.

4. Commit your inputs: `study.yaml`, `conditions/`, `analyses/` and
   `figures/`.

   ```bash
   git -C lipase_363K add -A
   git -C lipase_363K commit -m "Analyses for the paper"
   ```

5. Set a git user name and email, if git has none:

   ```bash
   git config --global user.name "Jane Doe"
   git config --global user.email "jane.doe@example.org"
   ```

## 3. Freeze

```bash
polyzymd study freeze lipase_363K
```

```
froze /home/me/lipase_363K as study-v1 (dc94243466cf)
manifest: 21 study files, 2 conditions, 10 replicates hashed
deposit: /home/me/lipase_363K/deposit; files to upload in /home/me/lipase_363K/deposit/upload
warning: metadata.doi is not set: reserve a DOI for the study in Zenodo, add it here and refreeze (deposit/UPLOAD.md says how)
next: follow /home/me/lipase_363K/deposit/UPLOAD.md, which says how to reserve the DOI, upload and publish on Zenodo; PolyzyMD uploads nothing
```

Freeze stops before it writes anything in these cases:

- An input file is not committed. The message names the files. Commit them
  and freeze again. The tag and the deposit hold only committed files, and
  the manifest must describe exactly those files.
- Git has no user name and email.

Otherwise, freeze does these steps:

1. It checks whether the stored results of each analysis still match the
   study. It compares the config hashes, the equilibration window, the
   stride, the code of each function, the settings and the PolyzyMD version.
   A stale run gives a warning that names what changed. Run it again with
   `polyzymd analyze RUN --study lipase_363K`.
2. It warns in these cases too:
   - The production lengths of the conditions differ by more than 10 %. A
     difference can then come from the simulated time. Analyze with `until`
     (see {ref}`study-until`), or make the short simulations longer.
   - The `progress.json` of a replicate lists production segments that are
     not on disk.
   - The config of a condition does not agree with its topology. Examples: a
     substrate residue that is missing from the topology, or polymer residues
     in a condition whose config has no polymers, or the reverse.
   - An OpenMM replicate records no OpenMM version. For a GROMACS study,
     freeze always asks you to state the GROMACS version in the methods.
   - A file that `build_manifest.json` lists is missing from the replicate
     folder.
3. It hashes every trajectory and topology file on this machine (SHA-256,
   computed once per file and cached).
4. It writes gzipped copies of the engine inputs of each replicate to
   `deposit/`: the OpenMM system XML and topology, or the GROMACS `.tpr`,
   `.top`, `.itp` and `.mdp`. It also writes the final frame of each
   replicate.
5. It writes these files to the study folder:

   | File | Holds |
   |---|---|
   | `manifest.json` | Each study file, trajectory, engine input and final frame, by size and SHA-256. The package versions. The full config of each condition. The production length of each replicate. The git commit. The warnings |
   | `CITATION.cff` | Citation File Format 1.2.0: the paper as `preferred-citation`, and PolyzyMD and the trajectory deposits under `references` |
   | `.zenodo.json` | Zenodo deposit metadata: `isSupplementTo` the paper, `requires` PolyzyMD, `references` the trajectories |
   | `md_checklist.yaml` | The reliability and reproducibility checklist of Communications Biology (2023), filled in from the manifest. Review each answer |
   | `system_summary.csv` | The box, atoms, waters, ions and composition of each replicate (checklist 4a) |

6. It commits these files and `results/`, and only those. It tags the commit
   `study-v1`, then `study-v2`, and so on. `--tag NAME` sets another tag.
7. It lays out `deposit/` and prepares the upload: `deposit/upload/`,
   `deposit/trajectories.csv` and `deposit/UPLOAD.md`.

More about the manifest and the checklist:

- The manifest records, for each replicate, the PolyzyMD and OpenMM versions
  that built and ran it (`simulated_with`, from `build_manifest.json` and
  `progress.json`).
- The manifest follows the JSON Schema `manifest-1.schema.json`. The schema
  ships with PolyzyMD, and freeze writes it into the deposit. `study.yaml`
  has its own schema, `study-1.schema.json`.
- Distance restraints in the config of a condition act in every phase. So
  item 3c of `md_checklist.yaml` reports that condition as restrained
  (biased) sampling. It lists the type, atoms, distance and force constant of
  each restraint.

## 4. Upload and publish with `deposit/UPLOAD.md`

PolyzyMD does not upload or publish. A Zenodo publication is permanent and
gets a DOI, so you do that step. `freeze` prepares the files in `deposit/`,
which git ignores:

| Path | Holds |
|---|---|
| `UPLOAD.md` | The steps for this study: reserve its DOI, add the files, fill in each Zenodo form field, review, publish, and make new versions |
| `upload/` | The files to add to the Zenodo upload, no more and no less. See below |
| `trajectories.csv` | Each trajectory and topology file by size and SHA-256, in batches that each fit one Zenodo record |
| `README.md` | Made from `metadata:`: what the study is and why, its authors, how to cite the paper, the dataset and PolyzyMD, its contents, how to reproduce it, and the verdict of each run. Your own `README.md` of the study stays in `study/` as you wrote it |
| `study/`, `engine_inputs/`, `final_frames/`, and the top-level files | The same content, not zipped, for inspection |

`upload/` holds these files:

- `README.md`, `CITATION.cff` and `manifest.json`, not zipped, so that Zenodo
  shows a preview;
- one zip file each for the study, the engine inputs and the final frames. A
  record holds at most 100 files, and Zenodo shows what is in a zip file.

No file in the deposit names a path of your machine. Records and reports name
files relative to the study, or to the folder of the replicate folders. The
deposited configs say `projects_directory: .` and `scratch_directory: data`.
A comment in them tells a reader how to point the study at a copy of the
trajectories. The config hash leaves out both directories, so stored results
still match.

Do these steps, which `UPLOAD.md` gives in full:

1. In Zenodo, select **New upload**. Answer that the upload has no DOI yet.
   Then select **Get a DOI now!**. Zenodo reserves a DOI for the draft. The
   DOI is registered when you publish, and lost if you delete the draft.
2. Set `metadata.doi` in `study.yaml` to that DOI. Commit.
3. Freeze again, so that the files contain their own DOI.
4. Add the files of the new `deposit/upload/` to the draft.
5. Fill in the form from the table in `UPLOAD.md`.
6. Review, and select **Publish**.

Zenodo fixes the files a short time after you publish. Its help pages say 30
and 45 days. For a later change, freeze again and upload a **New version**.
Each version has its own DOI. The concept DOI always points to the latest
version.

When a freeze has written `UPLOAD.md`, `polyzymd study check` prints
`publish: follow deposit/UPLOAD.md`.

## 5. Deposit the trajectories

PolyzyMD does not upload trajectories. `deposit/trajectories.csv` puts the
files of each condition into batches that fit one Zenodo record:

- 50 GB and 100 files by default;
- up to 200 GB with a quota increase from the storage settings of the draft.

Each replicate stays in one batch. The file names any single file over the
50 GB file limit of Zenodo.

1. Deposit each batch on Zenodo or in another data repository.
2. List each DOI under `metadata.related.trajectories`, with the conditions
   that it holds.
3. Commit, and freeze again.

## 6. Freeze again when the paper is published

1. Set `metadata.related.paper.doi`. Also set the trajectory DOIs, if you did
   not know them before.
2. Commit, and freeze again. The next tag has the new citations.
3. Upload the new deposit to Zenodo as a new version of the record. To add
   DOIs without changing the files, edit the metadata of the published
   record instead.

## Reproduce a published study

A reader downloads `study/` and the trajectories, and runs:

```bash
polyzymd study locate ~/Downloads/zenodo_1234567 --study lipase_363K --verify
```

```
No polymer: runs [1, 2, 3, 4, 5] under /home/me/Downloads/zenodo_1234567/no_polymer
No polymer: 10 files match manifest.json (SHA-256)
```

`--verify` checks the SHA-256 of each located file against `manifest.json`.
Without it, `locate` checks only the sizes. A changed or missing file is
named, and the command exits with code 2.

The reader can make the figures without the trajectories, from
`pz.Study("study.yaml").results(run)`. See {doc}`study_yaml`.

## Cite PolyzyMD

Each frozen study cites PolyzyMD in `CITATION.cff` and `.zenodo.json`. Its
generated README says how to cite PolyzyMD. Keep these citations when you
publish.
