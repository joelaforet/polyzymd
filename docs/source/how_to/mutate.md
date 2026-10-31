# Mutate a protein

```{index} !mutation, !point mutation, !non-canonical amino acid, !d-amino acid
```

Replace residues of a protein before it is built: with another of the 20
canonical amino acids, its D form, a non-canonical amino acid from the
Chemical Component Dictionary (CCD), or a side chain you write as SMILES.

## Step 1: Write the mutations

Mutations go in their own list under the protein, and PolyzyMD applies them
before any modification:

```yaml
    mutations:
      - {site: A:S1500, to: CYS}                       # a canonical residue
      - {site: A:L20, to: LEU, stereo: D}              # its D form
      - {site: A:S1381, to: {ccd_file: structures/4II.cif}}   # a CCD component
      - {site: A:K12, to: {smiles: "*CCCCN=[N+]=[N-]", resname: AZK}}   # a SMILES side chain
      - {chains: [A, C, E], residues: [T301], to: ASN} # the same change on three chains
```

- `to: CYS` takes the residue from the CCD files mBuild ships for the 20
  canonical amino acids.
- A non-canonical residue needs its CCD file. Download `<code>.cif` from the
  [wwPDB CCD](https://www.wwpdb.org/data/ccd) and point `ccd_file` to it.
- A SMILES side chain has `*` on the atom that bonds to CA.

A residue mutated here is named by its new code in later steps: after
`S1381 → 4II`, a click on that residue uses `site: A:4II1381`.

## Step 2: Mind protonation states

A CCD component arrives in the protonation state of its definition: Asp and
Glu as the neutral acids, His doubly protonated. The audit warns about these
(`modifications.protonation_state`). Write the step you want:

```yaml
    modifications:
      - deprotonate: {site: A:D45, atom_name: OD2}
```

## Step 3: Name the parameters of non-canonical residues

The 20 canonical residues and their D forms keep ff14SB. Any other residue is
a site: its side chain from CB outward takes the `noncanonical` parameter set,
and its backbone keeps ff14SB.

```yaml
    modification_parameters:
      noncanonical:
        force_field: ["openff-2.3.0.offxml"]
        charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

## Step 4: Check and build

```bash
polyzymd validate -c config.yaml --plans
polyzymd build -c config.yaml
```

The plan shows each mutation with the residue before and after. A mutation
written under `modifications:` is refused with a pointer to `mutations:`.
