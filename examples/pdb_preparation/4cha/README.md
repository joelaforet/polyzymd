# 4CHA OpenFF PDB ingestion proof of concept

This directory contains safe example scripts for investigating 4CHA
alpha-chymotrypsin PDB ingestion with OpenFF and PolyzyMD conventions.

The scripts are intentionally narrow. They are **not** a production-ready
universal PDB-preparation workflow.

## Files

| File | Purpose |
|---|---|
| `prepare_4cha.py` | Select one 4CHA enzyme copy, remove waters/heterogens, relabel protein chains to `A`, and add hydrogens only |
| `validate_openff.py` | Run structural checks and direct `Topology.from_pdb()` validation |
| `nterminal_cystine_substructure.json` | NCYX custom-substructure proof-of-concept mapping used by the example |

## Run the example

Download the raw PDB:

```bash
mkdir -p structures
curl -L https://files.rcsb.org/download/4CHA.pdb -o structures/4CHA.pdb
```

Prepare the narrow 4CHA example:

```bash
pixi run -e build python examples/pdb_preparation/4cha/prepare_4cha.py \
  structures/4CHA.pdb structures/4cha_chain_a.pdb
```

Validate directly with OpenFF:

```bash
pixi run -e build python examples/pdb_preparation/4cha/validate_openff.py \
  structures/4cha_chain_a.pdb
```

Test the private OpenFF custom-substructure proof of concept by preparing a copy
that renames the N-terminal cystine residue to `NCYX`, then validating with the
custom mapping:

```bash
pixi run -e build python examples/pdb_preparation/4cha/prepare_4cha.py \
  structures/4CHA.pdb structures/4cha_chain_a_ncyx.pdb \
  --rename-first-cys-to-ncyx
```

```bash
pixi run -e build python examples/pdb_preparation/4cha/validate_openff.py \
  structures/4cha_chain_a_ncyx.pdb \
  --custom-substructures examples/pdb_preparation/4cha/nterminal_cystine_substructure.json
```

## Caveats

- `Topology.from_pdb(..., _custom_substructures=...)` uses a private OpenFF API.
- The `NCYX` mapping is a targeted proof of concept for the 4CHA N-terminal
  cystine case and a candidate for upstream OpenFF discussion, not a stable
  public solution.
- Do not ignore charge mismatch errors. If validation fails, inspect and curate
  the structure instead of forcing PolyzyMD to continue.
- `polyzymd clean-pdb` does not guarantee OpenFF ingestion.

See the human troubleshooting guide at
`docs/source/how_to/troubleshoot_openff_pdb_ingestion.md` and the lookup
reference at `docs/source/reference/openff_pdb_ingestion.md`.
