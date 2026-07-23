# Fast mixed conjugation E2E fixture

Small committed inputs for the opt-in fast mixed checkpoints. The N-linked case
uses cleaned 1UBQ with an ABC 3-mer at `A:LYS63:NZ`, the original three-branch
SMILES glycan at `A:ASN60:ND2`, and one free SBMA 3-mer.

The O-linked cases use `glycam_G57321FI_CONECT.pdb`, an explicit-hydrogen,
GLYCAM-named monosaccharide copied without coordinate, atom-name, or connectivity
changes on 2026-07-21 from the repository maintainer's supplied validation corpus:
`Conjugation_Tests/structures/glycam_G57321FI_CONECT.pdb`. The `G57321FI`
identifier is retained as source provenance. Its reducing-end `ROH` oxygen and
hydrogen (serials 1 and 2) are removed and C1 (serial 4) is linked to Ser or Thr.
The maintainer supplied and authorized this test fixture for this repository;
this note makes no broader claim about upstream database licensing.
Its SHA-256 checksum is
`69ec434e6ae97f4356833917d411f7b3e4651dcb6b3e77b5e8d95d93dc3028b4`.

Generated artifacts are written only to pytest temporary directories.
