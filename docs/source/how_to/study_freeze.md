# Publish a study with `polyzymd study freeze`

Use this when a study's analyses are done and you want to deposit it with a
paper, so that anyone can redraw its figures from the deposit and rerun its
analyses from the trajectories. For creating the folder, see
{doc}`study_folder`; for what the deposit contains and why, see
{doc}`../explanation/study_folders`.

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
  doi: "10.5281/zenodo.NNNNNNN"            # the study's own DOI, reserved in Zenodo (step 3)
  title: "LipA with SBMA-EGMA copolymers at 363 K"
  description: "Lipase A simulated without polymer and with SBMA-EGMA copolymers."
  purpose: "To test whether the copolymer stabilises the lid at high temperature."
  keywords: [molecular dynamics, lipase, polymer, OpenMM]
  system_type: [protein, polymer]
  authors:
    - {family-names: Laforet, given-names: Joseph R., orcid: "0000-0000-0000-0000",
       affiliation: "University of Colorado Boulder"}
  license: {data: CC-BY-4.0, code: MIT}     # SPDX identifiers; match LICENSE-data and LICENSE-code
  funding:
    - {funder: "...", award: "...", funder_doi: "10.13039/..."}
  related:
    paper: {title: "...", status: in-preparation, doi: "10.XXXX/placeholder"}
    trajectories:                           # one entry per trajectory deposit
      - {doi: "10.5281/zenodo.1234567", conditions: [No polymer, SBMA 50%]}
    experimental:
      - {doi: "10.1021/...", description: "Measured activity at 363 K"}
  zenodo: {communities: [], access_right: open}
```

| Field | Why |
|---|---|
| `doi` | The study's persistent identifier (FAIR F1), written into `CITATION.cff`, `.zenodo.json` and the manifest; reserve it in Zenodo before publishing |
| `title`, `description`, `keywords`, `authors` | Findability (FAIR F2) and citation |
| `purpose` | The main purpose of the simulations, part of the minimum metadata of Amaro et al. (2025) |
| `system_type` | The kinds of molecule simulated, as the Communications Biology checklist asks (3a) |
| `license` | How others may reuse the data and code (FAIR R1.1) |
| `related.paper` | The paper's DOI, once it has one; until then `status: in-preparation` |
| `related.trajectories` | The DOIs of the trajectory deposits, which are uploaded on their own |
| `related.experimental` | Experiments the simulations are compared with (checklist 2a) |

An unknown key is refused with the nearest spelling. Anything missing is a
warning and a `TODO` placeholder, never a refusal: freeze early and refreeze
when the gaps are filled. `polyzymd study check` prints how many gaps remain.

## 2. Commit and freeze

Commit your inputs (`study.yaml`, `conditions/`, `analyses/`, `figures/`),
run every analysis (`polyzymd analyze --study study.yaml`), and, if the runs
have no recorded trajectory hashes, record them once with
`polyzymd hash-trajectories --study study.yaml` (freeze warns when they are
missing; see {doc}`study_folder`). Then:

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

Freezing:

1. Checks whether each analysis's stored results still match the study:
   config hashes, equilibration window, stride, the code of each function,
   the settings and the PolyzyMD version. A stale run is a warning naming
   what changed; rerun it with `polyzymd analyze RUN --study`. It also warns
   when:
   - the conditions' production lengths differ by more than 10%, so that a
     difference may come from simulated time (analyse with `until`, or
     extend the short runs);
   - a replicate's `progress.json` records production segments that are not
     on disk;
   - a condition's config disagrees with its topology: a substrate residue
     missing from the topology, or polymer residues in a condition whose
     config enables no polymers, or the reverse;
   - a replicate records no OpenMM version.
2. Hashes every trajectory and topology file on this machine (SHA-256,
   computed once per file and cached), and writes gzipped copies of each
   replicate's engine inputs (OpenMM system XML and topology, or GROMACS
   `.tpr`, `.top`, `.itp` and `.mdp`) and its final frame to `deposit/`.
3. Writes these files to the study folder:

   | File | Holds |
   |---|---|
   | `manifest.json` | Every study file, trajectory, engine input and final frame by size and SHA-256; package versions; each condition's fully resolved config; production length and frames analysed per replicate; the git commit; the warnings |
   | `CITATION.cff` | Citation File Format 1.2.0: the paper as `preferred-citation`; PolyzyMD and the trajectory deposits under `references` |
   | `.zenodo.json` | Zenodo deposit metadata: `isSupplementTo` the paper, `requires` PolyzyMD, `references` the trajectories |
   | `md_checklist.yaml` | The Communications Biology reliability and reproducibility checklist (2023), filled from the manifest; review each answer. Distance restraints in a condition's config act in every phase, so item 3c reports those conditions as restrained (biased) sampling, with each restraint's type, atoms, distance and force constant |
   | `system_summary.csv` | Box, atoms, waters, ions and composition of every replicate (checklist 4a) |

   The manifest also records, for each replicate, the PolyzyMD and OpenMM
   versions that built and ran it (`simulated_with`, from `build_manifest.json`
   and `progress.json`), and follows the JSON Schema `manifest-1.schema.json`,
   which ships with PolyzyMD and is written into the deposit; `study.yaml`
   has its own, `study-1.schema.json`.

4. Commits those files and `results/`, and only those, and tags the commit
   `study-v1` (then `study-v2`, ...; `--tag NAME` chooses). Freeze refuses to
   start while any input is uncommitted, and names those files: the tag and
   the deposit hold only committed files, and the manifest must describe
   exactly them. Commit your inputs first. Freeze also refuses when git has
   no user name and email.
5. Lays out `deposit/`, and prepares the upload: `deposit/upload/`,
   `deposit/trajectories.csv` and `deposit/UPLOAD.md`.

## 3. Upload and publish, by following `deposit/UPLOAD.md`

PolyzyMD uploads and publishes nothing: publishing on Zenodo is permanent and
mints a DOI, so that step is yours. `freeze` prepares everything for it in the
gitignored `deposit/`:

| Path | Holds |
|---|---|
| `UPLOAD.md` | The steps for this study: reserving its DOI, the files to add, the value of every Zenodo form field, reviewing and publishing, and new versions |
| `upload/` | Exactly the files to add to the Zenodo upload: `README.md`, `CITATION.cff` and `manifest.json` unzipped, so Zenodo previews them, and the study, the engine inputs and the final frames as one zip each, because a record holds at most 100 files and Zenodo shows what is inside a zip |
| `trajectories.csv` | Every trajectory and topology file by size and SHA-256, grouped into batches that each fit one Zenodo record |
| `README.md` | Written from `metadata:`: what the study is and why, its authors, how to cite the paper, the dataset and PolyzyMD, its contents, how to reproduce it, and each run's verdict. Your study's own `README.md` stays inside `study/` as written |
| `study/`, `engine_inputs/`, `final_frames/`, and the top-level files | The same content, unzipped, for inspection |

Nothing in the deposit names a path on your machine. Records and reports
name files relative to the study or the folder holding the runs, and the
deposited configs say `projects_directory: .` and `scratch_directory: data`
with a comment saying how a reproducer points the study at their copy of the
runs. The config hash does not include either directory, so stored results
still match.

The steps of `UPLOAD.md`, in short:

1. In Zenodo, choose **New upload**, answer that the upload has no DOI yet,
   and press **Get a DOI now!**. Zenodo reserves a DOI for the draft; it is
   registered when you publish and lost if you delete the draft.
2. Set `metadata.doi` in `study.yaml` to it, commit, and freeze again, so the
   files carry their own DOI.
3. Add the files of the new `deposit/upload/` to the draft, and fill in the
   form from the table in `UPLOAD.md`.
4. Review, and press **Publish**. Files are fixed shortly afterwards (Zenodo's
   help pages say 30 and 45 days); for a later change, freeze again and upload
   a **New version**. Each version has its own DOI, and the concept DOI always
   points to the latest.

`polyzymd study check` prints `publish: follow deposit/UPLOAD.md` once a freeze
has written it.

## 4. Trajectories

PolyzyMD does not upload trajectories either. `deposit/trajectories.csv` groups
each condition's run files into batches that fit one Zenodo record: 50 GB and
100 files by default, up to 200 GB with a quota increase from the draft's
storage settings. Replicates stay whole, and any single file over Zenodo's
50 GB file limit is named. Deposit each batch where it suits you, on Zenodo or
in a data repository; list each DOI under `metadata.related.trajectories`,
with the conditions it holds, and freeze again.

## 5. Refreeze when the paper is out

Set `metadata.related.paper.doi`, and the trajectory DOIs if they were not
known yet, commit, and freeze again: the next tag carries the updated
citations, and goes to Zenodo as a new version of the record (edit the
published record's metadata to add DOIs without changing files).

## Reproduce a published study

Someone reproducing the study downloads `study/` and the trajectories, and
runs:

```bash
polyzymd study locate ~/Downloads/zenodo_1234567 --study lipase_363K --verify
```

```
No polymer: runs [1, 2, 3, 4, 5] under /home/me/Downloads/zenodo_1234567/no_polymer
No polymer: 10 files match manifest.json (SHA-256)
```

`--verify` checks every located file's SHA-256 against `manifest.json`;
without it, only sizes are checked. A changed or missing file exits 2 and is
named. The figures can be redrawn without the trajectories, from
`pz.Study("study.yaml").results(run)` (see {doc}`study_yaml`).

## Cite PolyzyMD

Every frozen study cites PolyzyMD in `CITATION.cff` and `.zenodo.json`, and
its generated README says how to cite it. Please keep those citations when
you publish.
