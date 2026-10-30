# How modified residues are parameterized

```{index} !modified residue, !site, !capped model, !integer charge shift
```

A protein in PolyzyMD takes ff14SB, a force field fit to the 20 canonical
amino acids. A glycan, a dye, a grafted polymer or a non-canonical amino acid
has no ff14SB parameters. This page explains the rule PolyzyMD 1.4 uses for
those atoms, why, and what the audit checks.

## The rule

Each mutation to a non-canonical residue and each attachment defines a
**site**: the side chain of the changed residue from CB outward, plus every
atom of the attached group.

| Atoms | Bonded terms and Lennard-Jones | Charges |
|---|---|---|
| Residues without a change | ff14SB | ff14SB library charges |
| Backbone of a changed residue (N, H, CA, HA, C, O) | ff14SB | ff14SB library charges of that residue as it was before the change |
| The site | the parameter set named for its type (Sage 2.3.0, or GLYCAM06j for glycans) | NAGL on a capped model of the site, then a uniform shift to the integer charge; GLYCAM06j's own charges for glycans |

The boundary is the CA-CB bond. Phi, psi, the peptide bonds and the backbone
charges of a changed residue stay ff14SB, so its backbone is the same as its
neighbours'. Any term that touches a site atom (the CA-CB bond, the angles at
CA that involve CB, chi1, and everything further out) comes from the site's
force field.

**Charges of a site.** NAGL charges a small model of the site: the site's
residues plus two bonds of neighbours, each cut bond capped with a hydrogen.
The site's charges must then sum to the formal charge of the changed residue
and its group, minus the ff14SB charges of the backbone atoms that stay. A
uniform shift over the site atoms closes the gap. The build stops if the shift
exceeds `max_charge_shift` (0.005 e per atom by default), because a large
shift means the site is too small to absorb its charge without distorting it.

**Termini.** When a change touches a terminus (an N-terminal label, a
C-terminal amide), the untouched atoms of that residue take ff14SB's
main-chain template, as Amber caps a chain with ACE and NME. A terminus left
neutral without any group on it (NH2 or COOH) is refused: ff14SB has no
template for it, and protein termini are normally simulated charged at
physiological pH. Cap a truncated chain with ACE or NME instead; if you need
a neutral terminus (a low- or high-pH study), please
[open an issue](https://github.com/joelaforet/polyzymd/issues).

**Canonical changes.** A mutation to one of the 20 canonical residues, its D
form, or a protonation change to a state ff14SB defines (Lys to LYN, His to
HID, ...) makes no site. The residue stays fully ff14SB.

## Why the backbone stays ff14SB

A study on ubiquitin compared NAGL charges of a residue before and after three
modifications (an octanoyl and a methacrylate-like acyl on Lys63, an N-linked
GlcNAc on Asn60). The modification changed the backbone charges by at most
0.001 e. Replacing ff14SB's backbone charges with NAGL's, on the other hand,
changed them by up to 0.29 e (on CA) and 0.22 e on N from the difference
between the two charge models alone. Keeping ff14SB
on the backbone keeps a modified residue's backbone the same as its
neighbours'.

GLYCAM06j follows the same rule: its glycosylated Asn, Ser and Thr (NLN, OLS,
OLT) carry exactly the backbone charges of ff14SB's ASN, SER and THR.

## Why glycans use GLYCAM06j

A comparison on glycosylated ubiquitin (single sugars on Asn, Ser and Thr, and
chitobiose, a trisaccharide core and Man5 on Asn) put the Sage 2.3.0 + NAGL
rule against GLYCAM06j:

- **Charges were close.** The site charges differed from GLYCAM06j's by about
  0.08 e per atom on average, part of it because GLYCAM06j gives aliphatic
  hydrogens no charge (summed over each heavy atom and its hydrogens, 0.073 to
  0.087 e). Dipoles of single-sugar sites agreed within 0.4 D; sites with
  several sugars came out 2 to 2.2 D larger than GLYCAM06j's. NAGL and AM1-BCC
  differed by about 0.015 e per atom, and the integer shift stayed below
  0.002 e per atom.
- **Bonded terms were not.** Sage made the glycosidic and ring angles about
  2.4 times stiffer than GLYCAM06j, its torsion profiles differed by 6 to
  21 kJ/mol, and the minima of the Asn chi1 and of the Asn-glycan linkage
  torsion moved.

So glycans use GLYCAM06j. Sage with NAGL stays the general rule for the
groups no published force field covers, such as dyes, labels, grafted
polymers and non-canonical amino acids. A glycan on another force field
builds, with the warning `parameters.glycan_force_field`.

## GLYCAM06j's 1-4 pairs

ff14SB and Sage scale every 1-4 interaction by the same factors (1/1.2 for
electrostatics, 1/2 for Lennard-Jones). GLYCAM06j leaves selected 1-4 pairs of
a glycan unscaled: about 70 per sugar. Scaling them like the rest of the
protein changes the energy by about 150 kJ/mol for one sugar and about
1,080 kJ/mol for Man5.

OpenFF Interchange holds one 1-4 scale per system, so PolyzyMD builds the
glycan terms with OpenMM's own GLYCAM06j and writes each glycan 1-4 pair with
its own parameters: as a NonbondedForce exception for OpenMM, and as a
`[ pairs ]` line of function type 2 for GROMACS. Both engines match OpenMM's
native GLYCAM06j to about one part in a million. The workarounds this needs,
and the upstream issues behind them, are tracked in
[issue #281](https://github.com/joelaforet/polyzymd/issues/281).

## What the audit checks

The audit rebuilds each domain on its own and compares the built system with
it term by term. The checks specific to modified proteins:

| Check | Fails when |
|---|---|
| `modifications.declared` | the build's changes or parameter sets differ from the config's |
| `modifications.applied` | a planned bond, added or removed atom, or formal charge is missing from the built system |
| `modifications.sites` | a site holds atoms the rule does not put in it |
| `parameters.modified_backbone` | a changed residue's backbone is not ff14SB |
| `charges.site_integer` | a site's charges do not sum to its target |
| `charges.site_shift` | a site needed a shift over the bound |
| `parameters.one_four` | a 1-4 pair is not scaled as its owning force field scales it |

Every check is listed in {doc}`../reference/physics_checks`.

## Limits

- **Interim rule.** OpenFF's Rosemary is planned as one force field for
  proteins, small molecules and modified residues, which would remove the
  boundary. PolyzyMD offers its alpha as an opt-in ({doc}`../how_to/rosemary`).
- **Chi1.** With the boundary at CA-CB, the chi1 torsion of a changed residue
  comes from the site's force field, which was fit to small molecules. Where a
  side chain's orientation matters to a result, check chi1 against a QM scan
  of the capped residue.
- **Sage for other groups.** Sage's bonded terms have not been compared with
  a published force field for dyes, polymers or non-canonical residues.
