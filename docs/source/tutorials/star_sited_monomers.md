# Write star-sited monomer SMILES

Dynamic polymers are built from **star-sited SMILES**. Each monomer is given as
its repeat unit, as it is in the chain, with `*` on the two atoms that bond to
the neighbouring repeat units. Each chain end gets a cap with one `*`.

| What | Example | `*` sites |
|---|---|---|
| Methyl methacrylate (MMA), as sold | `C=C(C)C(=O)OC` | none |
| MMA repeat unit | `[*]CC([*])(C)C(=O)OC` | 2 |
| MMA as a chain end that kept its C=C (the chain ends PolyzyMD built before 1.3.0) | `[*]C=C(C)C(=O)OC` | 1 |
| Hydrogen cap | `[*][H]` | 1 |
| Methyl cap | `[*]C` | 1 |

The first `*` bonds toward the head of the chain, the second toward the tail.
For a vinyl, acrylate or methacrylate monomer, the two linking atoms are the
carbons of the C=C double bond, which is a single bond in the chain.

```{figure} ../_static/renders/star_sites.webp
:figclass: pz-render
:class: dark-light
:alt: Three panels of the same blue molecule: the SBMA monomer with its carbon-carbon double bond marked; the repeat unit with two small dark spheres for its two star sites; and a coiled chain of five units with halos along its backbone where the units are joined.

From SBMA to a chain, at one scale. Left: the monomer as sold, with the
two carbons of its C=C marked. Middle: the repeat unit
`[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]`, the C=C now single,
each `*` (dark) where the bond to the next unit points. Right: a 5-mer
that PolyzyMD built for `examples/polymer`, with the four bonds between
units marked; each one replaced the `*` of two units. Its end units keep
their C=C.
```

## The notebook

{download}`star_sited_monomers.ipynb <../../../src/polyzymd/templates/notebooks/star_sited_monomers.ipynb>` is a sandbox
for writing these SMILES. It:

1. draws the monomer SMILES you paste (SBMA by default), with every atom
   labelled by its index;
2. lists the C=C double bonds, the usual polymerizable bond;
3. puts `*` on the two atoms you name, redraws the repeat unit, and prints the
   star-sited SMILES;
4. writes the chain-end form of your monomer (the monomer with one `*` and
   its C=C kept, as at the ends of the chains PolyzyMD built before 1.3.0) and prints a
   `type: polymer` entry to paste into your config;
5. builds one short chain with the PolyzyMD builder (optional).

It does not change your config. Check each result by eye.

`polyzymd project init` puts a copy in your project's `notebooks/` folder.
Open it there in the PolyzyMD build environment (change `~/polyzymd` to the
folder of your PolyzyMD clone):

```bash
cd my_paper/notebooks
pixi run -e build --manifest-path ~/polyzymd/pixi.toml jupyter lab --no-browser star_sited_monomers.ipynb
```

JupyterLab prints a link with a token, such as
`http://localhost:8888/lab?token=...`. Open it in your web browser.

On Windows with WSL, `localhost` in the Windows browser may not reach WSL (the
page says the connection was refused). Then start JupyterLab with
`--ip=0.0.0.0` added, find the WSL address with `hostname -I`, and open
`http://<that address>:8888/lab?token=...` instead.

## Use the result

Copy the repeat unit into `smiles:` of each monomer and the caps, with their
probabilities, into `end_groups:`. {doc}`../how_to/dynamic_polymers` describes every key.
`polyzymd validate` reads each SMILES with RDKit and stops if a monomer does
not have exactly two `*` or a cap does not have exactly one.
