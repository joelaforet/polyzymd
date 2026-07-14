# Fast mixed conjugation E2E fixture

Small committed inputs for the opt-in fast mixed checkpoint. The runtime build
uses cleaned 1UBQ with an ABC 3-mer at `A:LYS63:NZ`, the validated three-branch
glycan at `A:ASN60:ND2`, and one free SBMA 3-mer. The glycan input is the full
three-branch SMILES checkpoint (`117` heavy atoms, `225` atoms with hydrogens,
formula `C64H108N2O51`), not a mono-NAG surrogate.

Generated artifacts are written only to pytest temporary directories.
