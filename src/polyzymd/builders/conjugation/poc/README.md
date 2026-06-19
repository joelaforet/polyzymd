# Conjugation POC Assets

Historical proof-of-concept assets for automatically conjugating NHS-ester
polymers to lysine residues on a protein. The old helper scripts that drove the
POC have been retired now that the public conjugation API is the maintained
entry point.

**GitHub issue**: <https://github.com/joelaforet/polyzymd/issues/51>

## Quick start

Use `public_api_conjugation_walkthrough.ipynb` for the maintained runnable
example. It exercises the supported public API rather than the retired POC
helper scripts.

`conjugation_poc_walkthrough.ipynb` and the bundled PDBs, images, cached
polymer files, and output snapshots are retained as historical reference
material only.

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
├── public_api_conjugation_walkthrough.ipynb ← maintained runnable example
├── conjugation_poc_walkthrough.ipynb        ← historical POC notebook
├── PEG_LYS_Schematic.png      ← historical schematic
├── *.pdb                      ← historical input and output structures
├── data/
│   └── ...                    ← retained POC data assets
└── output/                    ← retained generated artifacts and snapshots
```

The retired helper scripts were removed from this directory:

- `conjugation_poc.py`
- `generate_polymer.py`
- `free_polymer_placement.py`
- `pdb_convention.py`
- `run_conjugate_simulation.py`

## How it works

1. **Protein loading** — Load PDB with RDKit, extract ff14SB partial charges,
   identify reactive LYS-NZ sites.

2. **Polymer generation** — Generate an NHS-SBMA polymer with 3D coordinates
   through the maintained public API workflow.

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
