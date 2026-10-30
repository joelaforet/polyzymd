# Graft a polymer chain onto a protein

```{index} !graft, !polymer conjugate, !pegylation, !nhs ester
```

Bond a polymer chain to a protein residue, built from the same star-sited
monomers as the soluble chains of {doc}`dynamic_polymers`.

## Step 1: Define the monomers once

A grafted chain takes its monomers and caps from a `type: polymer` solute of
the same config (it must be dynamic). If the system has no soluble chains,
list `monomers` and `end_groups` in the fragment instead.

## Step 2: Choose the sequence

Positions count the whole chain, 1 to `length`, the two caps included. A
13-mer with an NHS-activated monomer in the middle and the other positions
drawn from the monomer probabilities is:

```yaml
length: 13
sequence: {fixed: {7: NHS}}
```

`sequence` also takes `mode: block` (runs such as `[{H: 1}, {SBMA: 5},
{NHS: 1}, {SBMA: 5}, {H: 1}]`, the caps included) or `mode: explicit` (every position listed). A monomer used only
at pinned positions may have probability 0.

## Step 3: Attach the chain

```yaml
    modifications:
      - deprotonate: {site: A:K8, atom_name: NZ}
      - attach:
          type: polymer
          site: A:K8
          atom_name: NZ
          # NHS ester aminolysis: the amine bonds the ester carbon, NHS leaves
          reaction: "[N:1]([H:2])[H].[C:3](=[O:4])[O][N]>>[N:1]([H:2])[C:3]=[O:4]"
          attach_at: 7
          fragment:
            polymer:
              monomers_from: pSBMA
              length: 13
              sequence: {fixed: {7: NHS}}
    modification_parameters:
      polymer:
        force_field: ["openff-2.3.0.offxml"]
        charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

`attach_at` is the position that bonds. For a chain grown from the protein or
a PEG with a reactive end, it is 1 (a cap position is allowed when that cap is
a monomer's end unit).

## Step 4: Decide how replicates differ

Each replicate draws its own chain, seeded by the replicate number and the
site, so a set of replicates samples the sequence space the way soluble
chains do. To give every replicate the same conjugate, set `random_seed` in
the `polymer` fragment. Each replicate's drawn sequence is in its build
record.

## What the build does

The whole grafted chain and the Lys side chain from CB outward form one site
(segment `P001`): Sage 2.3.0 parameters and NAGL charges on a capped model,
like any other attached group. Grafted chains count toward the polymer mass
fraction and are listed apart from soluble chains in
`composition.polymer_sequences`. To separate their contacts from those of the
soluble chains, see {doc}`analyze_modified_proteins`.
