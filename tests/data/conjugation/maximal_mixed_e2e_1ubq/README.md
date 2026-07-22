# Maximal mixed-conjugation E2E fixtures

Inputs unique to the dedicated opt-in maximal acceptance case. Generated build,
export, and MD artifacts are written only to the pytest temporary directory.

`glycam_G42666HT_CONECT.pdb` was copied byte-for-byte on 2026-07-22 from
`/home/joelaforet/Shirts-Lab-Linux/Conjugation_Tests/structures/glycam_G42666HT_CONECT.pdb`.
The source SHA-256 is
`09ac19da3795aad0078d79946d6b4fe39cf202797d73e3187d9ac824b12c4525`.
The only transformation was stripping trailing spaces from its 57 `ATOM`
records; coordinates, atom/residue names, record order, and `CONECT` records are
unchanged. The committed fixture SHA-256 is
`db28bbca394aee3edfa7fa6b1ac359377429942b22f3926b0d9122593084019a`.

The O-glycan reuses the repository's committed
`tests/data/conjugation/fast_mixed_e2e_1ubq/glycam_G57321FI_CONECT.pdb` fixture,
created on 2026-07-21 from
`/home/joelaforet/Shirts-Lab-Linux/Conjugation_Tests/structures/glycam_G57321FI_CONECT.pdb`.
The source SHA-256 is
`3cba7e8aae4ffd3060ad2f56f0900978e00ab20059edd8b3cc131e8c438eea91`.
The only transformation was stripping trailing spaces from its 30 `ATOM`
records; coordinates, atom/residue names, record order, and `CONECT` records are
unchanged. The committed fixture SHA-256 is
`69ec434e6ae97f4356833917d411f7b3e4651dcb6b3e77b5e8d95d93dc3028b4` and
the cleaned `1ubq.pdb` is reused from the same committed fast fixture.

The repository maintainer supplied the source corpus and explicitly authorized
these files for PolyzyMD validation. This provenance statement does not make a
broader claim about upstream database licensing.

The opt-in maximal-acceptance protocol is deliberately fixed at 293 K and 1 atm:
uncapped minimization; 20 ps restrained protein-backbone-heavy-atom NVT; 20 ps
unrestrained NPT; and 300 ps NPT production with a 2 fs LangevinMiddle timestep,
1/ps friction, MC barostat every 25 steps, CUDA mixed precision, DCD/state every
1,000 steps, checkpoint every 25,000 steps, and exact 150-frame readback. This
test protocol is not a universal recommendation for production MD.
