# Build Pipeline Optimization Summary

**Branch:** `patch/v1.0.1-bridges2`
**Date:** February 2026
**Test system:** Fibronectin 8-10 + SBMA/OEGMA copolymers (12-126 chains, 308K+ atoms)

---

## Problem Statement

The PolyzyMD system build pipeline (pre-simulation) was taking ~2.5 hours for a
308K-atom system (12 polymer chains, 7 unique types). Two bottlenecks dominated:

| Bottleneck | Duration | % of Total |
|-----------|----------|-----------|
| `Interchange.combine()` (serial left-fold) | 84 min | 55% |
| `to_openmm(combine_nonbonded_forces=False)` | 53 min | 35% |
| Everything else (Packmol, parameterization, PDB I/O) | 15 min | 10% |
| **Total** | **~152 min** | |

These timings are from SLURM logs on Bridges-2 HPC (RM partition, 128 cores,
512 GB RAM).

---

## Bottleneck #1: `Interchange.combine()` — 84 minutes

### Root Cause

The original code created one `Interchange` object per molecule instance (~1140
for a typical solvated system), then combined them one-by-one in a serial
left-fold. Each `combine()` call does `deepcopy` + topology merge + charge cache
rebuild — O(N^2) total cost.

An earlier optimization (commit `f5eb105`) improved this by batching by SMILES
type (one Interchange per unique type, ~7 batches), reducing combines from ~1140
to ~11. This brought parameterization from impossible to 84 minutes, but
`combine()` itself was still the dominant cost.

Tree-based (divide-and-conquer) combining was attempted but **does not work** —
`_DUPLICATE` parameter keys in OpenFF cause `UnsupportedCombinationError` when
combining two objects that both contain them (only left + right works, not
right + right).

### Solution: Single-call approach

**Key insight:** `ForceField.create_interchange()` already handles multi-molecule
topologies efficiently. Internally, `Topology.identical_molecule_groups` deduplicates
SMIRKS matching and charge assignment so work scales with the number of *unique*
molecule types (~12), not total molecules (~100K).

The fix passes the entire solvated topology to a single `create_interchange()`
call with pre-computed `charge_from_molecules` for all species. This eliminates
`combine()` entirely.

**Commit:** `a2af9af` — `perf: replace batched Interchange.combine() with single
create_interchange() call`

### Benchmark Result (local, 4 polymers, 199K atoms, box-padding 0.3)

| Step | Time |
|------|------|
| `create_interchange()` (single call) | 61.8s |
| `to_openmm_topology()` | 118.9s |
| `to_openmm(system)` | 176.1s |
| Total build | 356.8s (~6 min) |

For the full 308K-atom system on Bridges-2 HPC, `create_interchange()` is
projected at ~1-2 min (faster hardware, and the function scales sub-linearly
with molecule count due to deduplication).

---

## Bottleneck #2: `to_openmm()` — 53 minutes

### Root Cause (two independent issues)

**Issue A — Unnecessary `combine_nonbonded_forces=False`:**
Our code called `interchange.to_openmm(combine_nonbonded_forces=False)` since the
very first commit, with no documentation or reason. This flag triggers
`_create_multiple_nonbonded_forces()`, a slow code path that separates vdW, Coulomb,
and 1-4 interactions into individual OpenMM force objects. Our simulation runner
does nothing that requires split forces.

**Issue B — O(N^2) list membership in `_create_multiple_nonbonded_forces()`:**
Inside that slow code path, `openmm_pairs` was a Python `list` used for `in`
membership testing — O(N) per check. For 308K atoms (~340K exceptions x ~8K
1-4 pairs), this produced ~2.7 billion list scans.

### Solution

**Issue A:** Switch to `combine_nonbonded_forces=True` (the default). Single-line
change.

- Commit `5c0c198` — `fix: use combine_nonbonded_forces=True in to_openmm()`

**Issue B:** Upstream fix — convert `list` to `set` for O(1) lookups.

- Filed **openff-interchange#1444** (issue)
- Opened **openff-interchange#1445** (our PR)
- Matt Thompson merged the fix himself via **#1446**, released in **openff-interchange v0.5.1**
- Our monkey-patch (`openmm_patches.py`) was removed as dead code:
  - Commit `9aca3de` — removed monkey-patch usage from build pipeline
  - Commit `13eeefc` — deleted the module

---

## Bottleneck #3: PDB CONECT Record Overflow — Solvation Crash

### Root Cause

When the solute (protein + polymers) exceeds ~9,999 atoms, OpenMM's
`PDBFile.writeFile()` hex-encodes CONECT atom indices beyond 5 digits (the PDB
fixed-width format limit). Packmol's Fortran parser uses fixed-width reads and
cannot interpret hex-encoded indices, causing a crash during solvation.

This affects any system with ≥30 polymer chains (~10,485 solute atoms). The
126-chain target system has ~30,843 solute atoms with 26,643 CONECT records
that would overflow.

**Path:** `SolventBuilder.solvate()` → OpenFF `pack_box()` →
`_create_solute_pdb()` → OpenMM `PDBFile.writeFile()` → hex CONECT records →
Packmol crash.

### Solution

New `solvate_with_packmol()` function in `utils/packmol.py` replaces OpenFF's
`pack_box()` for the solvation step. Two-layer fix:

1. **Strip CONECT records** from the solute PDB after OpenMM writes it. Packmol
   only needs atomic coordinates; bond topology is preserved in the OpenFF
   `Topology` object.
2. **`ignore_conect` keyword** in Packmol input file for versions > 21.1.3 as
   belt-and-suspenders (our installed Packmol 21.1.3 doesn't support this, so
   CONECT stripping alone handles the fix).

**Commit:** `bd40354` — `fix: strip PDB CONECT records to prevent solvation
crash for large systems`

### Test Results

| Test | Solute Atoms | CONECT Stripped | Total Atoms | Solvation Time | Result |
|------|-------------|-----------------|-------------|----------------|--------|
| 4 chains | 5,077 | 877 | 199,664 | 95s | PASS |
| 30 chains | 10,485 | 6,285 | 220,425 | 101s | PASS |
| **126 chains** | **30,843** | **26,643** | **233,643** | **158s** | **PASS** |

The 30-chain test is the critical one — it crosses the ~10K atom threshold that
previously crashed. The 126-chain test validates the full target system.

---

## Fix #4: Polymer Shell Packing — Polymers Scattered Throughout Box

### Root Cause (two independent issues)

**Issue A — No shell constraint:** The Packmol input used `inside box` for
polymers with no exclusion zone, so polymers filled the entire box volume
uniformly instead of forming a coating shell around the protein.

**Issue B — Centering mismatch:** The protein's original (non-centered)
coordinates were used in the final output topology, while polymers were packed
around a *centered copy* of the protein. This spatial mismatch placed the
protein off to one side of the box with polymers clustered elsewhere.

### Solution

Added `outside box` constraint in Packmol input for every polymer structure
block, creating a rectangular exclusion zone equal to the protein's bounding
box ± Packmol tolerance. Combined with the existing `inside box` (outer
boundary), this creates a rectangular shell where polymers are packed.

The three-shell architecture:
- **Inner zone (protein bbox):** Exclusion zone — no polymers placed here
- **Shell zone (padding):** Polymers packed here, thickness = `packing.padding` (nm)
- **Outer zone:** Filled with water + ions by the solvation step

Fixed the centering bug by using `centered_solute` consistently in the output
topology. Forced `_use_pbc = False` in `pack_polymers()` since Packmol's PBC
mode does not support per-structure `outside box` constraints.

Additional improvements:
- `run_packmol()` accepts exit code 173 as non-fatal (logs warning, continues
  with best solution found)
- New `_max_molecule_diameter_angstrom()` helper computes largest polymer
  diameter for shell-thickness validation
- Emits a warning when shell thickness < polymer diameter, with a recommended
  `packing.padding` increase

**Commit:** `6db2ffc` — `fix: pack polymers in shell around protein and fix
centering mismatch`

### Key Parameters

- `packing.padding`: Shell thickness per side in nm. Default was 1.5 nm (too
  thin for ~28.5 Å diameter polymers). Changed to 3.0 nm in test config.
- The user noted that oversized boxes are acceptable because NPT equilibration
  with dynamic box resizing will shrink the box to the correct density.

### Test Results (simple exclusion zone, packing.padding=3.0 nm)

| Test | Packing Time | Solvation | Total Pipeline | Packmol Exit |
|------|-------------|-----------|----------------|-------------|
| 4 chains | 0.9s | 73.4s | 95.2s | 0 (success) |
| 30 chains | 3.5s | 132.3s | 165.2s | 0 (success) |
| **126 chains** | **49.5s** | **157.0s** | **236.1s (3.9 min)** | **0 (success)** |

All tests use `--box-padding 0.3` (minimal solvation) and `--skip-interchange`
(testing packing geometry only). Shell geometry for fibronectin 8-10:
- Protein bbox: ~43 × 71 × 105 Å
- Exclusion box: protein bbox ± 2.0 Å (Packmol tolerance)
- Shell thickness min: 25.5 Å (close to polymer diameter of 28.5 Å but viable)

---

## Fix #5: Chain ID Assignment Mismatch — Restraints on Wrong Atoms

### Root Cause

`_assign_pdb_identifiers()` in `system_builder.py` used a sequential `chain_idx`
counter that incremented only when a component was present. When substrate was
absent (`n_substrate_molecules == 0`), the substrate block was skipped entirely,
causing polymers to be assigned chain **B** and solvent to start at chain **C**.

But all downstream code hardcodes fixed chain conventions:
- `SystemComponentInfo`: `polymer_chain_id="C"`, `substrate_chain_id="B"`
- `AtomGroupResolver.from_topology()`: chain A → protein, B → substrate, C → polymer
- `_get_polymer_heavy()`: filters chain C for non-hydrogen elements

| Component | Actual PDB (buggy) | What code expects |
|-----------|-------------------|-------------------|
| Protein   | Chain A           | Chain A ✓         |
| Polymers  | Chain **B**       | Chain **C** ✗     |
| Solvent   | Chain **C**       | Chain **D** ✗     |

### Consequence

During equilibration, `polymer_heavy` restraints targeted **chain C atoms** —
which were actually water molecules. The log showed `Added positional restraints
to 9,999 atoms` for `polymer_heavy`. This is exactly 1/3 of 29,997 chain-C
water atoms (H₂O has 1 heavy atom per 3 total), not the expected 1,176 polymer
heavy atoms.

The 12-chain polymer system has:
- 2,556 total polymer atoms across 12 chains
- 1,176 heavy atoms (~46% heavy, matching SBMA/OEGMA SDF atom compositions)
- The previous 9,999 "polymer heavy atoms" was a dead giveaway of the bug

### Solution

Replace the sequential `chain_idx` counter with fixed chain letter assignments:

```python
PROTEIN_CHAIN = "A"
SUBSTRATE_CHAIN = "B"
POLYMER_CHAIN = "C"
SOLVENT_START_IDX = 3  # index of 'D'
```

Each component block uses its fixed letter directly. Absent components (e.g.,
no substrate) simply don't assign any atoms — but don't shift other components'
chain letters. This makes the `SystemComponentInfo` hardcoded defaults correct
by construction.

**Commit:** `16433b0` — `fix: enforce fixed chain ID assignment (A=protein,
B=substrate, C=polymer, D+=solvent)`

### Verification

12-chain fibronectin test (`--count 12 --box-padding 0.3`):

| Chain | Atom Count | Component |
|-------|-----------|-----------|
| A     | 4,200     | Protein ✓ |
| B     | 0 (reserved) | Substrate (absent) ✓ |
| C     | 2,556     | Polymers (12 chains) ✓ |
| D-J   | 202,854   | Solvent (water + ions) ✓ |

Polymer heavy atoms on chain C: **1,176** (correct, was 9,999 before fix).

---

## Projected Total Improvement

| Step | Before (Bridges-2) | After (projected Bridges-2) |
|------|--------|--------------------|
| Polymer packing (Packmol) | 3s | 3s |
| Solvation (Packmol) | 5 min 30s | 5 min 30s |
| PDB save | 46s | 46s |
| Parameterization (`create_interchange`) | 3 min 7s | ~1-2 min |
| Interchange combination (`combine()`) | **84 min 14s** | **0s** (eliminated) |
| OpenMM extraction (`to_openmm()`) | **53 min 4s** | **~2-3 min** |
| **Total** | **~152 min** | **~10-12 min** |

This is a **~13x speedup** on the build pipeline.

### Local Benchmark (4 chains, 199K atoms, box-padding 0.3)

Full end-to-end pipeline on local machine:

| Step | Time |
|------|------|
| Pack polymers | 1.0s |
| Solvation | 95.3s |
| Save PDB | 21.6s |
| `create_interchange()` | 65.7s |
| `to_openmm()` (topology + system + positions) | 288.2s |
| **Total** | **477.0s (~8 min)** |

Note: `to_openmm()` is 288s locally but would be ~3 min on Bridges-2 HPC with
faster CPUs. The local machine is ~3x slower for CPU-bound OpenFF operations.

---

## All Commits on `patch/v1.0.1-bridges2`

Listed newest-first, optimization-related commits only:

| Commit | Description |
|--------|-------------|
| `16433b0` | `fix:` Enforce fixed chain ID assignment (A=protein, B=substrate, C=polymer, D+=solvent) |
| `6db2ffc` | `fix:` Pack polymers in shell around protein and fix centering mismatch |
| `bd40354` | `fix:` Strip PDB CONECT records to prevent solvation crash for large systems |
| `a2af9af` | `perf:` Single-call `create_interchange()`, eliminate `combine()` entirely |
| `13eeefc` | `chore:` Remove unused `openmm_patches.py` module |
| `9aca3de` | `refactor:` Remove monkey-patch import/call from build pipeline |
| `5c0c198` | `fix:` Use `combine_nonbonded_forces=True` in `to_openmm()` |
| `33053b7` | `perf:` Monkey-patch `to_openmm()` list-to-set (now superseded by upstream fix) |
| `f5eb105` | `perf:` Cache-first polymer lookup, SMILES cache, Packmol sort, batch combining |

---

## Files Changed

| File | Changes |
|------|---------|
| `src/polyzymd/builders/system_builder.py` | New `_create_interchange_single_call()`, removed wrapper methods, `combine_nonbonded_forces=True`, LibraryCharges log suppression, `box_vectors_nm` parameter for shell packing, **fixed chain ID assignment** (sequential counter → fixed letters A/B/C/D+) |
| `src/polyzymd/builders/polymer.py` | Cache-first molecule lookup, SDF filename fix |
| `src/polyzymd/builders/solvent.py` | Replaced `pack_box()` with `solvate_with_packmol()` for CONECT overflow protection |
| `src/polyzymd/utils/packmol.py` | Descending-count sort, `solvate_with_packmol()`, `_strip_conect_records()`, `_check_ignore_conect_supported()`, shell packing via `inner_exclusion_box_angstrom` + `outside box` constraint, exit code 173 handling, `_max_molecule_diameter_angstrom()` |
| `src/polyzymd/utils/boxvectors.py` | New `get_topology_bbox_bounds()` helper for protein bounding box extraction |
| `src/polyzymd/utils/__init__.py` | Added `get_topology_bbox_bounds` export |
| `src/polyzymd/config/schema.py` | Added `box_vectors` field to `PolymerPackingConfig` |
| `src/polyzymd/utils/openmm_patches.py` | **Deleted** (monkey-patch no longer needed) |

---

## Upstream Contributions

| Item | Status |
|------|--------|
| openff-interchange#1444 (issue: `list` vs `set` in `to_openmm`) | Closed, fixed |
| openff-interchange#1445 (our PR) | Closed (author applied fix via #1446) |
| openff-interchange#1446 (Matt's PR, `list`->`set`) | Merged, released in v0.5.1 |

---

## Remaining Work

### Must Do: Full 126-chain Validation on Bridges-2 HPC (with Interchange)

The full build pipeline has been validated locally:
- **Shell packing (126 chains):** PASSED — polymers packed in shell, 49.5s, exit code 0
- **Solvation (126 chains):** PASSED — CONECT overflow fix works, 238K atoms, 2.6 min
- **Full pipeline (4 chains):** PASSED — solvation → Interchange → OpenMM, 477s
- **Chain ID assignment (12 chains):** PASSED — chain C=2,556 polymer atoms, 1,176 heavy

The 126-chain system with Interchange creation **should be tested on Bridges-2**
where the CPU-bound parameterization steps will be 2-3x faster. The profiling
harness (`profile_build.py`) in `BRIDGES_2_TESTING/fib_8_to_10/` is ready to use.

### Impact on Existing Bridges-2 Trajectory

The completed Bridges-2 trajectory at
`BRIDGES_2_TESTING/Fibronectin_8_to_10_SBMA-OEGMA_A75_B25_dynamic_test`
was run **before** the chain ID fix. During its equilibration:
- **Heating stage:** polymer heavy restraints were applied to water oxygens
  (chain C = water, not polymers). The polymer chains were unrestrained.
- **Polymer relaxation stage:** only protein_heavy was restrained (correct).
- **Free equilibration and production:** no restraints (correct).

The production trajectory may be scientifically usable despite the incorrect
equilibration restraints, but should be treated with caution and potentially
rerun with the fixed code.

### Nice to Have

- **Draft upstream feature request for custom Packmol constraints** — OpenFF's
  `pack_box()` does not support per-molecule Packmol constraints (`outside box`,
  `inside sphere`, etc.). Must draft with the user before filing.
- **File upstream issue for PDB CONECT overflow** — OpenFF's `pack_box` fails
  when solute exceeds 9999 atoms. Our workaround (CONECT stripping) is in
  `polyzymd`, but the root cause is in OpenFF's `_create_solute_pdb()` using
  OpenMM's `PDBFile.writeFile()` without CONECT suppression.
- **Investigate `to_openmm_topology()` slowness** — 102s for 200K atoms seems
  high. Likely an OpenFF internal issue (topology graph construction). Worth
  profiling on Bridges-2 where it may be faster.
- **Update openff-interchange to v0.5.1** in `polyzymd-bridges2` conda env to
  get the native `list->set` fix (currently we avoid the slow path entirely by
  using `combine_nonbonded_forces=True`, so this is not blocking).

---

## Upstream Issue Drafting Notes

These notes provide context for drafting issues/feature requests on upstream
repositories. **Framing is critical:** these are NOT bugs — the upstream code
works correctly. Frame everything as "performance improvement" or "enhancement"
using collaborative, respectful language. The upstream maintainers (especially
Matt Thompson at OpenFF) have been responsive and helpful.

### Our Info

- **GitHub user:** `joelaforet`
- **Our repo:** `github.com/joelaforet/polyzymd`
- **Branch with fixes:** `patch/v1.0.1-bridges2`
- **Upstream repo for Interchange:** `github.com/openforcefield/openff-interchange`
- **Already completed upstream:** #1444 (issue), #1445 (our PR, closed),
  #1446 (Matt's fix, merged in v0.5.1)

### Issue Draft 1: PDB CONECT Record Overflow in `pack_box()`

**Target repo:** `openforcefield/openff-interchange`

**Affected code path:**
`openff.interchange.components._packmol._create_solute_pdb()` calls
`openmm.app.PDBFile.writeFile()` which writes CONECT records. When the solute
exceeds ~9,999 atoms, OpenMM hex-encodes atom indices in CONECT records (PDB
format only allows 5-digit atom serial numbers). Packmol's Fortran parser reads
fixed-width columns and cannot parse hex-encoded indices, causing a crash.

**Reproduction:** Any call to `pack_box()` where the solute topology has >9,999
atoms. In our case: protein (4,200 atoms) + 30 polymer chains (6,285 atoms) =
10,485 atoms. The Packmol error is a Fortran read error on the CONECT lines.

**Our workaround:** We strip CONECT records from the PDB file after OpenMM
writes it, before Packmol reads it. Packmol only needs ATOM/HETATM records for
coordinate placement; bond topology is preserved in the OpenFF `Topology` object
and never needs to round-trip through PDB.

**Possible upstream fixes:**
1. Use `openmm.app.PDBFile.writeFile(..., keepIds=False)` — but this may not
   suppress CONECT. Needs investigation.
2. Strip CONECT records in `_create_solute_pdb()` after writing.
3. Use a format other than PDB for the solute (e.g., XYZ), though this may
   require Packmol changes.

### Issue Draft 2: Custom Packmol Constraints in `pack_box()`

**Target repo:** `openforcefield/openff-interchange`

**Current limitation:** `pack_box()` hard-codes `inside box 0. 0. 0. Lx Ly Lz`
for every molecule. There is no way to pass per-molecule Packmol constraints
such as `outside box` (exclusion zone), `inside sphere`, `outside sphere`,
`fixed`, or `nloop`. The function signature accepts `Molecule` objects and
counts but has no constraint parameter.

**Relevant code:** `openff.interchange.components._packmol.pack_box()` builds
Packmol input via string concatenation. Each molecule block gets only
`inside box`. The `build_packmol_input()` internal function is not exposed.

**Use case:** Enzyme-polymer conjugate simulations require polymers packed in
a **shell around the protein**, not scattered throughout the box. This needs:
- `inside box` (outer boundary — the full packing box)
- `outside box xmin ymin zmin xmax ymax zmax` (inner exclusion — protein bbox)

Packmol's `outside box` semantics: each atom must satisfy
`x < xmin OR x > xmax OR y < ymin OR y > ymax OR z < zmin OR z > zmax`.
Combined with `inside box`, this creates a rectangular shell.

**Our workaround:** We wrote our own `build_packmol_input()` and
`pack_polymers()` in `polyzymd/utils/packmol.py` that accept an
`inner_exclusion_box_angstrom` parameter. This works but means we maintain
a parallel Packmol wrapper.

**Proposed API sketch (minimal change):**
```python
# Option A: Per-molecule constraint strings
pack_box(
    molecules=[protein, polymer_a, polymer_b],
    number_of_copies=[1, 10, 5],
    # New parameter: list of extra Packmol constraint lines per molecule
    extra_constraints=[
        ["fixed 0. 0. 0. 0. 0. 0."],  # protein is fixed
        ["outside box 28.4 25.5 26.8 75.2 100.7 136.2"],  # polymer_a
        ["outside box 28.4 25.5 26.8 75.2 100.7 136.2"],  # polymer_b
    ],
)

# Option B: Global constraint dict
pack_box(
    molecules=[...],
    packmol_keywords={"nloop": 200},
    per_molecule_constraints=[...],
)
```

The key principle: expose Packmol's constraint system without reinventing it.
Users who need advanced packing (shells, spheres, fixed molecules) can pass
raw Packmol constraint strings.

---

## Polymer Cache Symlinks

Four symlinks were created in the test data's `.polymer_cache/` to bridge naming
convention mismatches between generated sequences and cached SDF files:

```
SBMA_seq=AAAAA         -> SBMA-OEGMA_seq=AAAAA
OEGMA-SBMA_seq=BAAAB   -> SBMA-OEGMA_seq=BAAAB
OEGMA-SBMA_seq=BAABB   -> SBMA-OEGMA_seq=BAABB
OEGMA-SBMA_seq=BABBB   -> SBMA-OEGMA_seq=BABBB
```

These are only needed for >12 chain counts where additional unique sequences are
generated beyond the base 7 types.

---

## How to Reproduce Benchmarks

```bash
cd /home/joelaforet/Desktop/enzyme_immobilization/BRIDGES_2_TESTING/fib_8_to_10/

# Quick local test (4 chains, small box)
mamba run -n polyzymd-bridges2 python profile_build.py \
    --config testing_key_config.yaml \
    --count 4 \
    --box-padding 0.3 \
    --output-dir _benchmark_quick

# A/B comparison (single-call vs batched)
mamba run -n polyzymd-bridges2 python benchmark_interchange.py \
    --config testing_key_config.yaml \
    --count 4 \
    --box-padding 0.3 \
    --skip-batched  # or --skip-single for the other path
```
