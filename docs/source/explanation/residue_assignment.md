# Chains and residue numbers in a built system

`polyzymd build` gives every atom of the system a chain ID and a residue
number. The chain ID tells you what the molecule is. The residue number
tells you which unit of that molecule the atom belongs to. Analyses and
restraints select atoms with these two identifiers, so you need to know how
the builder sets them.

## One residue is one repeat unit

The builder makes each residue one unit that you would select or count:

| Component | One residue is |
|---|---|
| Protein | One amino acid |
| Substrate or ligand | The whole molecule |
| Polymer | One monomer |
| Water, ion or co-solvent | One molecule |

Some topology tools give all molecules of one type the same residue
identity. Then you cannot tell two waters apart, or count the contacts of one
monomer. The builder sets the identifiers when it still knows which atoms
form which molecule. They cannot be rebuilt reliably from a finished
trajectory.

## Chain IDs

The chain letters are fixed. A component that is absent leaves its letter
unused, so the letters of the other components never move.

| Chain | Contents |
|---|---|
| A | Every protein molecule |
| B | The substrate or ligand |
| C | Every polymer chain |
| D, E, ... | Water, ions and co-solvents |

Therefore `chainid A` always selects the protein, and `chainid C` always
selects all polymer. The PDB format allows 9999 residues in a chain. The
builder puts solvent molecules 1 to 9999 on chain D, the next 9999 on chain
E, and so on. If the solvent needs more chains than the alphabet has after
D, the build stops with an error. Solvent chain IDs never repeat A, B or C.

## Protein residue numbers

The builder renumbers the protein residues from 1, in the order of the input
PDB. It does not keep the input numbers. If the first residue of a crystal
structure is number 5, it becomes residue 1, and every other residue moves
down by 4.

A protein with several molecules, such as a homodimer, is all on chain A. The
numbering continues from one molecule to the next. If each monomer has 100
residues, the first monomer has residues 1 to 100 and the second has 101 to
200.

Use the built numbers in restraints, `core:` and `regions:` selections and
residue lists. To compare with the input numbering, add the offset back. See
{doc}`../how_to/restraints`.

## Polymer residue numbers

Each monomer of each polymer chain is one residue on chain C. The numbering
does not restart for each chain. It continues from the last monomer of one
chain to the first monomer of the next. With two chains of 10 monomers, the
first chain has residues 1 to 10, and the second has 11 to 20.

Therefore a residue number on chain C identifies one monomer of one chain.
Residue 1 is the first monomer of the first chain only. To select the first
monomer of every chain, use the chain lengths to find the numbers. The
residue name is the monomer's name. A selection such as
`chainid C and resname <name>` selects every monomer of that type.

## Substrate and solvent residue numbers

The substrate is residue 1 of chain B. Water, ions and co-solvent molecules
are numbered from 1 on each solvent chain, one number per molecule. A solvent
molecule is identified by its chain and its residue number together.

## What the identifiers make possible

Because each unit has its own identifier, an analysis can answer questions
such as these:

- Which monomers touch the protein most often?
- Which waters stay near the active site?
- Which co-solvent molecules bridge the protein and the ligand?

The `contacts` analysis, for example, reports a value for each monomer type
by its residue name on chain C. A viewer can show the unit that an analysis
names, from its chain and residue number.

## For contributors

Code that adds a component or changes the builder must keep these rules:

- Keep the fixed chain letters. Do not shift them when a component is absent.
- Keep one residue per repeat unit. Do not merge molecules into one residue.
- Refer to a unit by its chain and residue number together.
