# Hydrogen bonds implementation verification

This page records how `polyzymd analyze hydrogen_bonds` chooses donors and
acceptors, why it does not use MDAnalysis's charge-based guess, and what its
values were checked against.

## What the code does

1. The loader reads the partial charges and covalent bonds of the OpenMM
   system that PolyzyMD saves beside each segment's trajectory,
   `<segment>_system.xml`: the charges of the `NonbondedForce`, the bonds of
   the `HarmonicBondForce`, and the constraints, which OpenMM uses for bonds
   to hydrogen and rigid water, without the hydrogen-hydrogen constraint of
   water. The topology PDB alone records bonds only for non-standard
   residues: on a lipase with an SBMA-EGMA polymer, 7,720 bonds for 78,410
   atoms, all in the polymer.
2. `functions.hbond_atoms` chooses the donatable hydrogens, those bonded to
   N, O or S, and the acceptors, every O and every N or S with at most two
   bonded atoms.
3. `functions.hydrogen_bonds` gives these atoms to MDAnalysis
   `HydrogenBondAnalysis` as selections, which pairs each hydrogen with its
   donor through the bonds and finds, on every frame, the donor-acceptor
   pairs within `d_a_cutoff` with the minimum image of the box and a
   donor-hydrogen-acceptor angle of at least `d_h_a_angle_cutoff`. Bonds
   within one residue are left out.

## Why not MDAnalysis's charge-based guess

`HydrogenBondAnalysis` can guess donors and acceptors from partial charges:
a donor is an atom bonded to a hydrogen with a charge below -0.5, and an
acceptor any atom with a charge below -0.5. On the same lipase system, with
the charges of its force field, that guess:

| | Atoms |
|---|---|
| Misses as donors | every backbone amide N-H (charge -0.42), Lys NZ, His ND1-HD1, Trp NE1 |
| Takes as acceptors | the quaternary ammonium nitrogen of SBMA (-0.60), Arg NE, NH1 and NH2, Asn ND2, some backbone N |
| Misses as acceptors | the ether oxygens of EGMA (-0.42), Met SD |

It found 69 donors, all in Arg, Asn, Gln, Ser, Thr and Tyr side chains. A
nitrogen with four bonds, or three bonds in an amide or guanidinium group,
has no lone pair to accept a hydrogen bond, whatever its charge, so
PolyzyMD chooses by element and valency instead. Its donors are those of
MDTraj's Baker-Hubbard code, N-H and O-H, with S-H added; its acceptors are
restricted to atoms that keep a lone pair, where MDTraj takes every N and O. On the lipase system
the rule gives 318 donatable hydrogens and 1,314 acceptors among the protein
and polymer atoms, including His NE2 (unprotonated), Met SD, and every
carbonyl, carboxylate, hydroxyl, ether, ester and sulfonate oxygen. The
chosen atoms are recorded in every report under `provenance.settings`.

## Agreement with the plugin used before this version

The plugin used before this version took every N and O of the groups as a
donor and acceptor candidate, paired hydrogens with donors within 1.2 Å, and
used 3.0 Å. On 716 production frames of a 363 K *B. subtilis* lipase A
replicate with a 50:50 SBMA-EGMA polymer, in September 2026:

| Calculation | Protein-polymer bonds per frame | Residue pairs per frame |
|---|---|---|
| The plugin | 8.849162011173185 | 8.19413407821229 |
| `hydrogen_bonds` with the plugin's selections and 3.0 Å | 8.849162011173185 | 8.19413407821229 |
| `hydrogen_bonds` with the element and valency rule and 3.0 Å | 8.849162011173185 | |
| `hydrogen_bonds` with the element and valency rule and 3.5 Å (default) | 13.406424581005586 | 12.018156424581006 |

The two selections give the same protein-polymer bonds at 3.0 Å on this
replicate: the acceptors the plugin added, such as the ammonium nitrogen,
never met the geometry. The default cutoff of 3.5 Å counts about half again
as many bonds.

## What the tests check

`tests/analyses/test_hydrogen_bonds_functions.py` checks, on small synthetic
systems, the distance and angle limits, both directions between two groups,
bonds within one residue, the element and valency rule, the counts against an
independent calculation, lifetimes by residue and atom pairs, occupancies,
pair labels, and the results, settings and refusals of `polyzymd analyze
hydrogen_bonds`; `tests/analyses/test_force_field_enrichment.py` checks the
reading of the system XML.

## Scope

The rule needs the bonds of every hydrogen of the groups; without the system
XML, or a topology with bonds, the analysis refuses to run unless the
hydrogens and acceptors are given as selections. The geometry is
MDAnalysis's, with the donor-acceptor distance and not the hydrogen-acceptor
distance. The counts depend on the cutoffs, which should be reported with the
results.

## See Also

- [Hydrogen Bonds Quick Start](../how_to/hydrogen_bonds.md) — commands and settings
- [How long polymer contacts last](analysis_contact_lifetimes.md) — the lifetime estimator
- [Shipped analysis functions](../reference/analysis_functions.md) — what each function measures
