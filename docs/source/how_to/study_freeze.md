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

## 1. Fill in the metadata

Add a `metadata:` block to `study.yaml`. These fields make the deposit
findable and citable:

```yaml
metadata:
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
run every analysis (`polyzymd analyze --study study.yaml`), then:

```bash
polyzymd study freeze lipase_363K --zip
```

```
froze /home/me/lipase_363K as study-v1 (dc94243466cf)
manifest: 21 study files, 2 conditions, 10 replicates hashed
deposit: /home/me/lipase_363K/deposit; zip /home/me/lipase_363K/deposit/lipase_363K-study-v1.zip
warning: the paper DOI is missing or a placeholder; refreeze once it is known
next: upload the trajectories and the files in deposit/ (see ...)
```

Freezing:

1. Checks whether each analysis's stored results still match the study:
   config hashes, equilibration window, stride, the code of each function,
   the settings and the PolyzyMD version. A stale run is a warning naming
   what changed; rerun it with `polyzymd analyze RUN --study`.
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

4. Commits those files and `results/`, and only those, and tags the commit
   `study-v1` (then `study-v2`, ...; `--tag NAME` chooses). Your uncommitted
   inputs are never committed: they are listed as a warning and are not part
   of the tagged study.
5. Lays out `deposit/` for upload.

## 3. Upload

`deposit/` is gitignored, and holds what to upload, file by file:

| Path | Holds |
|---|---|
| `manifest.json`, `CITATION.cff`, `.zenodo.json`, `README.md` | At the top, so the deposit is indexed and citable |
| `study/` | The tagged study, as `git archive` gives it |
| `engine_inputs/<condition>/replicate_<n>/` | Gzipped engine inputs |
| `final_frames/<condition>/` | Gzipped final-frame PDBs |
| `<study>-<tag>.zip` | With `--zip`, all of the above in one file |

Upload the files themselves rather than only the zip: archives whose
contents cannot be indexed hide data from search (MDverse, Tiemann et al.
2024). Upload the trajectories as their own deposits and put their DOIs in
`metadata.related.trajectories`. Zenodo reads `.zenodo.json` when a release
comes from GitHub; for a deposit made by hand, copy its fields into the
upload form.

## 4. Refreeze when the paper is out

Set `metadata.related.paper.doi`, and the trajectory DOIs if they were not
known yet, commit, and freeze again: the next tag (`study-v2`) carries the
updated citations, and Zenodo can take it as a new version of the deposit.

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
