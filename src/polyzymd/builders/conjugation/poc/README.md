# Conjugation POC — Polymer-Protein Covalent Linkage

Proof-of-concept for automatically conjugating NHS-ester polymers to lysine
residues on a protein, building a mixed force-field system (ff14SB + OpenFF
Sage), and running a full MD simulation workflow.

**GitHub issue**: <https://github.com/joelaforet/polyzymd/issues/51>

## Quick start

```bash
# 1. Install the conjugation-poc pixi environment (from repo root)
pixi install -e conjugation-poc

# 2. Open the Marimo notebook (interactive)
pixi run -e conjugation-poc conjugation-poc

# 3. Or run the full simulation workflow (headless)
pixi run -e conjugation-poc conjugation-sim
```

> **Note**: The simulation script (`run_conjugate_simulation.py`) imports the
> notebook (`conjugation_poc.py`) as a module to build the conjugate, then
> proceeds through solvation, equilibration, and production MD. It takes
> ~90 min on CPU.

## Test system

| Property | Value |
|---|---|
| Protein | Lipase A (*Bacillus subtilis*), 181 residues, 2721 atoms |
| Conjugation sites | LYS 23, LYS 44 (selected from eligible: 23, 35, 44) |
| Polymer | SBMA–NHS–SBMA (sequence `BAB`; A=NHS, B=SBMA) |
| Chemistry | NHS ester + LYS-NZ → amide bond + NHS leaving group |
| Protein FF | ff14SB (`ff14sb_off_impropers_0.0.4.offxml`) |
| Polymer FF | OpenFF Sage 2.0 (`openff-2.0.0.offxml`) |
| Water model | TIP3P |
| Ionic strength | 137 mM NaCl |

## File layout

```
poc/
├── README.md                  ← this file
├── __init__.py                ← package docstring
├── conjugation_poc.py         ← Marimo notebook: build conjugate + parameterize
├── generate_polymer.py        ← subprocess helper: polymer generation via polymerist
├── run_conjugate_simulation.py← full simulation: solvate → equilibrate → production
├── data/
│   └── NH3_terminal_His_proton_updated.pdb   ← Lipase A protein structure
└── output/                    ← simulation artifacts (git-ignored)
    └── .gitkeep
```

## How it works

1. **Protein loading** — Load PDB with RDKit, extract ff14SB partial charges,
   identify reactive LYS-NZ sites.

2. **Polymer generation** — Call `generate_polymer.py` as a subprocess (uses
   the `build` pixi env with mbuild/polymerist) to create an NHS-SBMA polymer
   with 3D coordinates.

3. **Placement** — Rigid-body translate/rotate each polymer to position the
   NHS carbonyl carbon near the target LYS-NZ, respecting bond angle geometry
   and avoiding steric clashes with the protein.

4. **Graph surgery** — RDKit `CombineMols` + `AddBond` to form the amide
   linkage. Remove NHS leaving group and two NZ hydrogens. Protein coordinates
   are kept immutable.

5. **Charge assignment** — NAGL charges on capped product-state polymer;
   ff14SB charges on protein atoms. Charge redistribution to maintain net
   integer charge.

6. **Parameterization** — `ForceField([ff14SB, Sage]).create_interchange()`
   with the combined OpenFF `Topology` and merged charge vector.

7. **Minimization** — Local energy minimization of the conjugate in vacuum
   (linkage neighborhood only, then full system).

8. **Solvation** — `solvate_topology()` with TIP3P water and 137 mM NaCl in
   a cubic box.

9. **Equilibration** — Three stages:
   - NVT heating 60→310 K (0.3 ns), protein+polymer heavy atoms restrained
   - NVT polymer relaxation at 310 K (0.5 ns), protein heavy atoms restrained
   - NPT free equilibration at 310 K (0.5 ns)

10. **Production** — 1 ns NPT at 310 K, 1 atm.

## Environment variables (optional overrides)

All paths default to locations relative to this directory. Override with:

| Variable | Default | Description |
|---|---|---|
| `CONJUGATION_PROTEIN_PDB` | `data/NH3_terminal_His_proton_updated.pdb` | Protein PDB input |
| `CONJUGATION_ASSEMBLED_PDB` | `output/conjugate_assembled.pdb` | Pre-minimization output |
| `CONJUGATION_MINIMIZED_PDB` | `output/conjugate_minimized.pdb` | Post-minimization output |
| `CONJUGATION_BUILD_PYTHON` | `<repo>/.pixi/envs/build/bin/python` | Python for polymer generation subprocess |
| `CONJUGATION_GENERATE_SCRIPT` | `generate_polymer.py` | Polymer generation script |
| `CONJUGATION_POLYMER_CACHE` | `.polymer_cache/` | Cached polymer structures |

## Validated results

From the completed test run (CPU, ~85 min):

- **47,508 atoms**: 2721 protein + polymer + 14,848 waters + 75 ions
- **Box**: 7.84 nm cube
- **Production**: T = 310.2 ± 1.3 K, density = 0.9995 ± 0.0023 g/mL
- All energies finite and stable throughout

## Known issues / limitations

- PDB loading loses bond orders → requires `DetermineBondOrders()` for NAGL
- RDKit vs OpenFF formal charge mismatch on ARG/HIS (~0.95e redistributed)
- `CombineMols` preserves stale conformers → must `RemoveAllConformers()`
  before `Molecule.from_rdkit()` and re-add the correct conformer
- CUDA platform validation requires a test `Context` (not just
  `getPlatformByName`); this machine falls back to CPU
