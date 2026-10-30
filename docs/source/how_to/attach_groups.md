# Attach a dye, a label or another group

```{index} !attach, !conjugation, !click chemistry, !thiol-maleimide, !dye, !label
```

Bond a small molecule to a protein residue: a fluorescent dye, a biotin, a
lipid tail or any group you can write as SMILES.

## Step 1: Write the group

Give the group as SMILES with `*` on the atom that bonds to the protein, and a
residue name of up to three characters:

```yaml
fragment: {smiles: "*C(=O)CCCCCCC", resname: OC8}   # an octanoyl group
```

A group can also come from a PDB file with `CONECT` records
(`{pdb_path: structures/group.pdb}`).

## Step 2: Prepare the residue

PolyzyMD changes nothing you did not write. A Lys must lose a proton before
its NZ can form an amide, so write that step first:

```yaml
    modifications:
      - deprotonate: {site: A:K63, atom_name: NZ}
      - attach:
          type: generic
          site: A:K63
          atom_name: NZ
          fragment: {smiles: "*C(=O)CCCCCCC", resname: OC8}
    modification_parameters:
      generic:
        force_field: ["openff-2.3.0.offxml"]
        charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

`examples/modified_proteins/config_octanoyl.yaml` is this config.

## Step 3: Pick how the bond forms

Without `reaction`, the `*` atom bonds to `atom_name` and one hydrogen leaves
each side. With a reaction, the reaction decides which atoms bond, which
leave and which bonds change order. mBuild names these reactions:

| `reaction` | Protein atom | Group |
|---|---|---|
| `amide coupling` | Lys NZ (deprotonated first) | a free carboxylic acid |
| `thiol-maleimide` | Cys SG | a maleimide |
| `azide-alkyne triazole` | the azide N3 of a residue such as AzF (`4II`) | an alkyne, such as DBCO |
| `N-glycosylation`, `O-glycosylation` | Asn ND2, Ser OG, Thr OG1 | a glycan (see {doc}`glycosylate`) |

Any other chemistry is a reaction SMARTS with the protein side first. An NHS
ester reacting with a Lys amine (NHS leaves) is:

```yaml
reaction: "[N:1]([H:2])[H].[C:3](=[O:4])[O][N]>>[N:1]([H:2])[C:3]=[O:4]"
```

A click that closes a ring forms two bonds. Add `merge: true` to make the
group part of the residue instead of a residue of its own:

```yaml
      - attach:
          type: generic
          site: A:4II1381
          atom_name: N3
          reaction: azide-alkyne triazole
          merge: true
          fragment: {smiles: "CC(=O)N1Cc2ccccc2C#Cc2ccccc21", resname: DBC}
```

## Step 4: Check and build

```bash
polyzymd validate -c config.yaml --plans
polyzymd build -c config.yaml
```

The plan shows the atoms that bond and leave and the change in formal charge.

## How the group is parameterized

The group and the side chain of its residue, from CB outward, form one site.
The site takes the force field you name in `modification_parameters`
(Sage 2.3.0 above) and NAGL charges computed on a small capped model of the
site. The backbone of the residue keeps ff14SB, like every other residue. A
small uniform shift brings the site to its integer charge; a shift over
`max_charge_shift` (0.005 e per atom) stops the build.
{doc}`../explanation/modified_residues` explains the rule.

If one entry names a parameter set for `generic` and the residue is a
non-canonical residue with its own `noncanonical` set (a click onto AzF), the
two sets must be the same.
