# OpenFF PDB Ingestion Rules

Use these rules when troubleshooting protein or PDB ingestion failures that reach
`openff.toolkit.Topology.from_pdb()` through PolyzyMD's enzyme builder.

## Living error log

When an agent or subagent encounters and diagnoses a new OpenFF PDB ingestion
error signature, the task is not complete until the durable context is updated or
the user explicitly defers the documentation update.

Update both:

- `docs/source/how_to/troubleshoot_openff_pdb_ingestion.md`
- `docs/source/reference/openff_pdb_ingestion.md`

Record:

1. the exact error text or shortest unique traceback excerpt
2. the likely cause
3. a diagnostic snippet or command
4. the acceptable fix
5. caveats, including when a workaround is private, experimental, or structure
   specific

Do not silently leave newly diagnosed failures in chat history, scratch notes, or
temporary scripts only.

## Required triage pattern

1. Reproduce the failure with direct OpenFF validation before changing PolyzyMD:

   ```python
   from openff.toolkit import Topology

   Topology.from_pdb("prepared_enzyme.pdb")
   ```

2. Inspect the input PDB as chemistry, not just as text. Check chains, TER
   records, residue names, disulfide connectivity, missing-coordinate records,
   heavy atoms, hydrogens, and terminal residues.
3. Treat OpenFF charge-mismatch diagnostics as real chemistry/connectivity
   problems. Do not ignore them and do not weaken validation just to continue.
4. Fix the structure upstream or provide a narrow documented proof of concept.
   Do not monkeypatch OpenFF internals, PolyzyMD builders, or toolkit residue
   templates in production code.

## Non-negotiable rules

- **Do not ignore charge mismatch.** A mismatch usually means the perceived
  residue graph, hydrogens, termini, formal charges, or bonds differ from the
  template OpenFF matched.
- **Do not monkeypatch.** Avoid patching OpenFF private functions or replacing
  toolkit templates at runtime. If a private API is used for a proof of concept,
  document it as private and experimental.
- **Do not claim `polyzymd clean-pdb` guarantees OpenFF ingestion.** It can help
  with simple cleanup, but OpenFF chemistry assignment is a separate validation
  step.
- **Keep imports lazy in scripts.** Example and diagnostic scripts should import
  OpenFF, OpenMM, PDBFixer, or MDAnalysis inside functions so basic module import
  does not require the full simulation stack.

## Disulfides, CYX, and connectivity

- OpenFF may treat `CYX` as a cysteine alias during PDB parsing, but `CYX` is not
  a stable public residue-template replacement for all disulfide cases.
- A disulfide cysteine SG atom should not carry an HG proton.
- `SSBOND` records alone are not always enough for every parser path. Verify
  SG-SG connectivity; curated PDBs commonly include appropriate `CONECT` records.
- If a cysteine is terminal, especially N-terminal, verify both terminal
  hydrogens and disulfide state. This is a common source of charge mismatch.

## Hydrogens, termini, and missing residues

- Hydrogens must be explicit before OpenFF PDB ingestion succeeds for proteins.
- Terminal residues need chemically consistent terminal atom names and hydrogen
  counts. Fix naming/protonation upstream rather than deleting atoms to silence
  errors.
- Missing-coordinate records (`REMARK 465`) and missing heavy atoms are modeling
  decisions. Do not let PolyzyMD invent missing residues automatically unless the
  user has explicitly approved the modeling choice and reviewed the result.

## Custom substructures

PolyzyMD exposes `enzyme.custom_substructures_path`, which is passed to
`Topology.from_pdb(..., _custom_substructures=...)`. The leading underscore means
this OpenFF argument is private/experimental. Use it only as a documented proof
of concept or targeted bridge while preparing an upstream OpenFF issue or pull
request.

The JSON shape is:

```json
{
  "RESNAME": {
    "[SMARTS:1]": ["ATOM1"]
  }
}
```

## 4CHA N-terminal cystine case

The raw 4CHA alpha-chymotrypsin workflow can pass simple structural checks after
selecting one biological copy, relabeling protein atoms to chain `A`, preserving
TER records, removing waters, and adding hydrogens, while still failing direct
OpenFF validation around the N-terminal disulfide cysteine and nearby residues.
This is an OpenFF chemistry-ingestion problem, not proof that PolyzyMD config is
wrong.

Use `examples/pdb_preparation/4cha/` as a proof-of-concept reference. The NCYX
custom-substructure example is structure specific and private-API based; treat it
as candidate evidence for an upstream OpenFF fix, not a universal production
preparation workflow.
