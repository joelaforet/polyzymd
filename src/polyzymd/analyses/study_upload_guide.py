"""Prepare a frozen study for upload to Zenodo, and say how to upload it.

PolyzyMD does not upload or publish anything: publishing is permanent and
mints a DOI, so it stays the author's step. ``polyzymd study freeze`` calls
:func:`prepare_upload`, which writes into ``deposit/``:

- ``upload/``: exactly the files to add to a Zenodo upload, within Zenodo's
  limits: the top-level metadata files, which stay previewable, and the
  study tree, engine inputs and final frames as one zip each (Zenodo
  recommends zipping more than 20 files and can browse inside a zip);
- ``trajectories.csv``: every trajectory and topology file the manifest
  lists, with its size and SHA-256, grouped into batches that each fit one
  Zenodo record;
- ``UPLOAD.md``: the steps for this study, from reserving its DOI to
  publishing, with the value of every Zenodo form field.

The Zenodo limits below are those documented on help.zenodo.org ("Manage
files", "Manage storage quota") and set in Zenodo's ``invenio.cfg``, as
checked in October 2026.
"""

from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

#: Bytes one Zenodo record holds by default.
RECORD_BYTES = 50_000_000_000
#: Bytes one record can hold after a self-service quota increase.
RECORD_BYTES_INCREASED = 200_000_000_000
#: Files one Zenodo record can hold; a quota increase does not change it.
RECORD_FILES = 100
#: Largest single file Zenodo accepts.
FILE_BYTES = 50_000_000_000
UPLOAD = "upload"
GUIDE = "UPLOAD.md"
TRAJECTORIES = "trajectories.csv"
#: Top-level files of the upload, kept out of the zips so they stay previewable.
TOP_FILES = ("README.md", "CITATION.cff", "manifest.json", "manifest-1.schema.json")


@dataclass(frozen=True)
class Batch:
    """Trajectory files of one or more replicates of a condition that fit one Zenodo record."""

    number: int
    condition: str
    replicates: tuple[str, ...]
    files: tuple[dict[str, Any], ...]

    @property
    def size(self) -> int:
        return sum(int(f["size"]) for f in self.files)


def _gb(size: float) -> str:
    """Return ``size`` bytes in decimal units, as Zenodo states its limits."""
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if size >= scale:
            return f"{size / scale:.2f} {unit}"
    return f"{int(size)} B"


def trajectory_batches(manifest: dict[str, Any]) -> tuple[list[Batch], list[str]]:
    """Group the manifest's trajectory files into batches that each fit one Zenodo record.

    A replicate's files stay together, and replicates of one condition are
    added to a batch until the next would pass :data:`RECORD_BYTES` or
    :data:`RECORD_FILES`. A replicate too large for one record is a batch of
    its own, and every file larger than :data:`FILE_BYTES` is named in the
    returned notes.
    """
    batches: list[Batch] = []
    notes: list[str] = []
    for label, condition in manifest.get("conditions", {}).items():
        current: list[tuple[str, list[dict[str, Any]]]] = []
        for index, replicate in sorted(
            condition.get("replicates", {}).items(), key=lambda kv: int(kv[0])
        ):
            files = [dict(f, condition=label, replicate=index) for f in replicate.get("files", [])]
            for item in files:
                if int(item["size"]) > FILE_BYTES:
                    notes.append(
                        f"{label} replicate {index}: {item['path']} is {_gb(int(item['size']))}, "
                        f"larger than Zenodo's {_gb(FILE_BYTES)} file limit; split it or deposit it "
                        "elsewhere"
                    )
            size = sum(int(f["size"]) for _, group in current for f in group)
            count = sum(len(group) for _, group in current)
            added = sum(int(f["size"]) for f in files)
            if current and (size + added > RECORD_BYTES or count + len(files) > RECORD_FILES):
                batches.append(_batch(len(batches) + 1, label, current))
                current = []
            current.append((index, files))
        if current:
            batches.append(_batch(len(batches) + 1, label, current))
    for batch in batches:
        if batch.size > RECORD_BYTES:
            notes.append(
                f"batch {batch.number} ({batch.condition} replicate {', '.join(batch.replicates)}) is "
                f"{_gb(batch.size)}: more than one record's default {_gb(RECORD_BYTES)}; increase the "
                f"record's quota (up to {_gb(RECORD_BYTES_INCREASED)}) or deposit it elsewhere"
            )
    return batches, notes


def _batch(number: int, label: str, group: list[tuple[str, list[dict[str, Any]]]]) -> Batch:
    return Batch(
        number, label, tuple(i for i, _ in group), tuple(f for _, files in group for f in files)
    )


def _zip(source: Path, target: Path) -> Path | None:
    if not source.is_dir() or not any(source.rglob("*")):
        return None
    return Path(
        shutil.make_archive(
            str(target.with_suffix("")), "zip", root_dir=source.parent, base_dir=source.name
        )
    )


def upload_folder(deposit: Path, name: str) -> list[Path]:
    """Write ``deposit/upload/``: the top-level files and one zip each of the study, engine inputs and final frames."""
    folder = deposit / UPLOAD
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    files = []
    for top in TOP_FILES:
        if (deposit / top).is_file():
            files.append(Path(shutil.copy2(deposit / top, folder / top)))
    for part, label in (
        ("study", f"{name}.zip"),
        ("engine_inputs", "engine_inputs.zip"),
        ("final_frames", "final_frames.zip"),
    ):
        made = _zip(deposit / part, folder / label)
        if made:
            files.append(made)
    return files


def _form_rows(
    zenodo: dict[str, Any], doi: str | None, licenses: dict[str, str] | None = None
) -> list[tuple[str, str]]:
    """Return the Zenodo upload form's fields with the values from ``.zenodo.json``."""
    creators = "; ".join(
        c["name"]
        + (f" (ORCID {c['orcid']})" if c.get("orcid") else "")
        + (f", {c['affiliation']}" if c.get("affiliation") else "")
        for c in zenodo.get("creators", [])
    )
    related = "; ".join(
        f"{r['relation']} {r['identifier']}"
        + (f" ({r['resource_type']})" if r.get("resource_type") else "")
        for r in zenodo.get("related_identifiers", [])
    )
    rows = [
        (
            "Digital Object Identifier",
            f"Yes, I already have one: {doi}" if doi else "No: press Get a DOI now (step 1)",
        ),
        ("Resource type", "Dataset"),
        ("Title", zenodo.get("title", "")),
        ("Publication date", zenodo.get("publication_date", "")),
        ("Creators", creators),
        ("Description", zenodo.get("description", "").replace("\n", " ")),
        (
            "Licenses",
            f"{licenses['data']} (data, results and figures) and {licenses['code']} (code in "
            "analyses/ and figures/); add both"
            if licenses
            else zenodo.get("license", "").upper(),
        ),
        ("Keywords and subjects", ", ".join(zenodo.get("keywords", []))),
        ("Version", zenodo.get("version", "")),
        ("Related works", related or "none"),
        ("Method (additional description)", zenodo.get("method", "")),
    ]
    if zenodo.get("communities"):
        rows.append(("Communities", ", ".join(c["identifier"] for c in zenodo["communities"])))
    if zenodo.get("grants"):
        rows.append(("Funding", ", ".join(g["id"] for g in zenodo["grants"])))
    return rows


def write_guide(
    deposit: Path,
    *,
    study_name: str,
    tag: str | None,
    meta: dict[str, Any],
    zenodo: dict[str, Any],
    upload: list[Path],
    batches: list[Batch],
    notes: list[str],
    warnings: list[str],
) -> Path:
    """Write ``deposit/UPLOAD.md``, the steps to publish this frozen study on Zenodo."""
    from polyzymd.citation import citation_line

    doi = meta.get("doi")
    files = "\n".join(f"| `upload/{p.name}` | {_gb(p.stat().st_size)} |" for p in upload)
    form = "\n".join(
        f"| {field} | {value} |" for field, value in _form_rows(zenodo, doi, meta.get("license"))
    )
    batch_rows = (
        "\n".join(
            f"| {b.number} | {b.condition} | {', '.join(b.replicates)} | {len(b.files)} | {_gb(b.size)} |"
            for b in batches
        )
        or "| - | no trajectories were on this machine when the study was frozen | | | |"
    )
    gaps = "\n".join(f"- {w}" for w in warnings) or "- none"
    trajectory_notes = "\n".join(f"- {n}" for n in notes)
    step1 = (
        f"The study's DOI is already set: `{doi}`. Use **Yes, I already have one** and enter it, "
        "or open the draft that reserved it."
        if doi
        else """1. On https://zenodo.org, choose **New upload**.
2. Under *Digital Object Identifier*, answer **No** and press **Get a DOI now!**. Zenodo
   reserves a DOI for this draft; it is registered only when you publish, and is
   lost if you delete the draft.
3. Put the DOI in `study.yaml`:

   ```yaml
   metadata:
     doi: "10.5281/zenodo.NNNNNNN"
   ```

4. Commit, and run `polyzymd study freeze` again: the next tag carries the DOI in
   `CITATION.cff`, `.zenodo.json` and `manifest.json`. Then upload that freeze's
   `deposit/upload/` to the same draft."""
    )
    text = f"""# Publish {study_name} on Zenodo

Written by `polyzymd study freeze` for `{tag or "an untagged freeze"}`. PolyzyMD prepares
the files and this guide; uploading and publishing are yours, because publishing
on Zenodo is permanent and mints a DOI. Test on https://sandbox.zenodo.org first
if you like: it is a separate account, and nothing there is real.

## Gaps to fill first

These are the warnings of the freeze. None stops you, but each one is missing
from the deposit:

{gaps}

## 1. Reserve the study's DOI

{step1}

## 2. Add the files

Add every file in `deposit/upload/` to the upload, no more and no less:

| File | Size |
|---|---|
{files}

The README, citation and manifest stay unzipped so Zenodo previews them; the
study tree, engine inputs and final frames are zipped because a record holds at
most {RECORD_FILES} files, and Zenodo shows what is inside a zip.

## 3. Fill in the form

Every value comes from `.zenodo.json`, which `freeze` wrote from `metadata:` in
`study.yaml`:

| Zenodo field | Value |
|---|---|
{form}

## 4. Review and publish

Preview the record, check the files and the citation, then press **Publish**.
After publishing:

- the DOI is registered and permanent, and the record can be deleted only within
  30 days, leaving a tombstone;
- files can no longer be changed after a short period (Zenodo's help pages give 30
  and 45 days); metadata can be edited at any time;
- to change files later, edit `study.yaml`, freeze again and upload the new
  `deposit/upload/` as a **New version** of the record. Each version gets its own
  DOI; the record's concept DOI always resolves to the latest. Cite the version
  DOI for reproducibility.

## Trajectories

PolyzyMD does not upload trajectories. `deposit/{TRAJECTORIES}` lists every trajectory
and topology file the manifest records, with its size and SHA-256, grouped into
batches that each fit one Zenodo record ({_gb(RECORD_BYTES)} and {RECORD_FILES} files by default; up
to {_gb(RECORD_BYTES_INCREASED)} with a quota increase from the draft's storage settings):

| Batch | Condition | Replicates | Files | Size |
|---|---|---|---|---|
{batch_rows}

{trajectory_notes}

Deposit each batch as its own record, on Zenodo or in a repository that suits
data this size, keeping the paths of `{TRAJECTORIES}` under one folder per condition.
Then list each DOI under `metadata.related.trajectories` in `study.yaml`,
with the conditions it holds, and freeze again; edit the study record's metadata
(Related works) to add them. Anyone who downloads them runs
`polyzymd study locate DOWNLOAD_DIR --verify`, which checks every file against
`manifest.json`.

## Automating it

For scripted uploads, Zenodo's REST API is documented at
https://developers.zenodo.org and https://inveniordm.docs.cern.ch; the files in
`deposit/upload/` and the fields of `.zenodo.json` are what it needs.

## Cite PolyzyMD

`CITATION.cff` and `.zenodo.json` already cite PolyzyMD; please keep them:

> {citation_line()}
"""
    path = deposit / GUIDE
    path.write_text(text)
    return path


def write_trajectories(deposit: Path, batches: list[Batch]) -> Path:
    """Write ``deposit/trajectories.csv``: batch, condition, replicate, path, size and SHA-256 of every file."""
    path = deposit / TRAJECTORIES
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["batch", "condition", "replicate", "path", "size_bytes", "sha256"])
        for batch in batches:
            for item in batch.files:
                writer.writerow(
                    [
                        batch.number,
                        item["condition"],
                        item["replicate"],
                        item["path"],
                        item["size"],
                        item["sha256"],
                    ]
                )
    return path


def prepare_upload(
    deposit: Path,
    *,
    study_name: str,
    tag: str | None,
    manifest: dict[str, Any],
    zenodo: dict[str, Any],
    warnings: list[str],
) -> dict[str, Path]:
    """Write ``deposit/upload/``, ``deposit/trajectories.csv`` and ``deposit/UPLOAD.md``; return their paths."""
    upload = upload_folder(deposit, f"{study_name}-{tag or 'untagged'}")
    batches, notes = trajectory_batches(manifest)
    trajectories = write_trajectories(deposit, batches)
    guide = write_guide(
        deposit,
        study_name=study_name,
        tag=tag,
        meta=manifest.get("metadata", {}),
        zenodo=zenodo,
        upload=upload,
        batches=batches,
        notes=notes,
        warnings=warnings,
    )
    return {"upload": deposit / UPLOAD, "trajectories": trajectories, "guide": guide}


def report_summary(path: Path) -> str:
    """Return a stored report's verdicts on one line, saying first when the report is partial."""
    import json

    try:
        report = json.loads(path.read_text())
    except (OSError, ValueError):
        return "no stored report"
    text = " ".join(report.get("verdict", []))
    if report.get("status", "complete") != "complete":
        problems = "; ".join(report.get("problems") or [])
        text = f"PARTIAL REPORT ({problems}). " + text
    return text


def deposit_readme(
    *,
    study_name: str,
    tag: str | None,
    meta: dict[str, Any],
    analyses: dict[str, Any],
    root: Path,
    project: bool = False,
) -> str:
    """Return the README at the top of the deposit, written from the study's metadata.

    It says what the study is and why it was run, who made it, how to cite it
    and PolyzyMD, what the deposit holds, how to reproduce it at each level,
    and what each analysis run found, as the verdict lines of its stored
    report.
    """
    import json

    from polyzymd.citation import citation_line

    authors = "; ".join(
        p.get("name")
        or ", ".join(x for x in (p.get("family-names"), p.get("given-names")) if x)
        + (f" (ORCID {p['orcid']})" if p.get("orcid") else "")
        for p in meta.get("authors", [])
    )
    paper = meta.get("related", {}).get("paper", {})
    doi = meta.get("doi")
    lines = [
        f"# {meta.get('title') or study_name}",
        "",
        str(meta.get("description", "")),
        "",
        f"**Purpose:** {meta.get('purpose', '')}",
        "",
        f"**Authors:** {authors}",
        "",
        f"**Version:** {tag or 'untagged'}" + (f" · **DOI:** {doi}" if doi else ""),
        "",
        "## How to cite",
        "",
        "Cite the paper, this dataset, and PolyzyMD, which produced the analyses "
        "(`CITATION.cff` holds all three):",
        "",
    ]
    if paper.get("title") or paper.get("doi"):
        lines.append(
            f"- Paper: {paper.get('title', '')}"
            + (f", doi:{paper['doi']}" if paper.get("doi") else "")
        )
    lines.append(f"- Dataset: {meta.get('title') or study_name}" + (f", doi:{doi}" if doi else ""))
    lines.append(f"- PolyzyMD: {citation_line()}")
    lines += [
        "",
        "## Contents",
        "",
        "| File | Holds |",
        "|---|---|",
        "| `manifest.json` | Every file by size and SHA-256, software versions, each condition's resolved config, and the production analysed |",
        "| `CITATION.cff` | How to cite the paper, this dataset and PolyzyMD |",
        "| `manifest-1.schema.json` | The JSON Schema `manifest.json` follows |",
        (
            f"| `{study_name}-{tag or 'untagged'}.zip` | The project: `project.yaml`, shared "
            "`analyses/`, `stats/` and `figures/`, and one folder per study (`study.yaml`, "
            "`conditions/`, `structures/`, `results/`) |"
            if project
            else f"| `{study_name}-{tag or 'untagged'}.zip` | The study: `study.yaml`, "
            "`conditions/`, `analyses/`, `figures/`, `results/` |"
        ),
        "| `engine_inputs.zip` | Each replicate's serialized engine inputs |",
        "| `final_frames.zip` | Each replicate's final frame |",
        "",
        "Trajectories are deposited separately; `manifest.json` lists each file with its SHA-256.",
        "",
        "## Reproduce",
        "",
        "Install the PolyzyMD version recorded in `manifest.json` (`versions`), unzip the "
        + ("project" if project else "study")
        + ", then:",
        "",
        "1. **Figures, without trajectories:** run the scripts in `figures/`, which read",
        '   `pz.Project(".").results(run)`.'
        if project
        else '   `pz.Study("study.yaml").results(run)`.',
        "2. **Analyses, from the trajectories:** download them, run",
        "   `polyzymd study locate DOWNLOAD_DIR --verify --study <study>` for each study, then",
        "   `polyzymd analyze --project .`."
        if project
        else "   `polyzymd study locate DOWNLOAD_DIR --verify`, then `polyzymd analyze --study study.yaml`.",
        "3. **Simulations:** build and run each `conditions/<name>/config.yaml` with PolyzyMD;",
        "   the replicate number is the random seed, so results agree within MD noise.",
        "",
    ]
    if project:
        # The project's results are each study's; the Studies section that
        # follows gives them.
        lines += [
            f"Licences: {meta.get('license', {}).get('data', 'CC-BY-4.0')} for data, results "
            f"and figures; {meta.get('license', {}).get('code', 'MIT')} for code.",
            "",
        ]
        return "\n".join(lines)
    lines += ["## Results", ""]
    for run in analyses:
        lines.append(f"- **{run}:** " + report_summary(root / "results" / run / "report.json"))
    lines += [
        "",
        f"Licences: {meta.get('license', {}).get('data', 'CC-BY-4.0')} for data, results and figures; "
        f"{meta.get('license', {}).get('code', 'MIT')} for code.",
        "",
    ]
    return "\n".join(lines)
