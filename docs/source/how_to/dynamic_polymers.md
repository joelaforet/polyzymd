# Generate polymers from SMILES

Build polymer chains from monomer SMILES with `generation_mode: "dynamic"`.
In this mode you do not need an SDF file for each sequence. PolyzyMD builds
each chain with Polymerist and ATRP (atom-transfer radical polymerization)
reaction templates, and assigns partial charges with the method in `charger`.

For the keys that both modes share, such as `monomers`, `length`, `count` and
`packing`, see {doc}`polymers`.

:::{admonition} Environment Setup
:class: tip

Run the build commands on this page in the `build` environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

## Choose a mode

| Mode | Use it when |
|---|---|
| `cached` | You have checked SDF files for each sequence and want to use exactly those structures. |
| `dynamic` | You have only the monomer SMILES, or you test a new monomer. |

## Write the `polymers:` section

```yaml
polymers:
  enabled: true
  generation_mode: "dynamic"
  type_prefix: "SBMA-EGPMA"
  reactions:
    initiation: "default"
    polymerization: "default"
    termination: "default"
  monomers:
    - label: "A"
      probability: 0.7
      name: "SBMA"
      residue_name: "SBM"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
    - label: "B"
      probability: 0.3
      name: "EGPMA"
      residue_name: "EGM"
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
  length: 5
  count: 2
  charger: "nagl"
  max_retries: 10
  cache_directory: ".polymer_cache"
```

`dynamic` mode adds these keys and rules:

| Key | Meaning |
|---|---|
| `smiles` | The SMILES of the monomer before polymerization, with its C=C double bond. Every monomer needs one. |
| `name` | The monomer name. Dynamic mode uses it to name the fragments, such as `SBMA_1-site` and `SBMA_2-site`, so give one for every monomer. `polyzymd validate` does not check this; the build stops instead. |
| `residue_name` | The 3-character residue name of the monomer in the topology. Optional. If you leave it out, PolyzyMD makes one from `name`. |
| `reactions` | The three ATRP reaction templates. `"default"` selects the templates that ship with PolyzyMD. You can also give the path of your own `.rxn` file. |
| `charger` | The partial-charge method for the chains: `nagl` (default), `am1bcc` or `espaloma`. |
| `max_retries` | The number of attempts to build a chain without a ring piercing (default 10). |
| `cache_directory` | The folder for the generated fragments and chains (default `.polymer_cache`). |
| `length` | Must be 3 or more. Polymerist needs at least one middle monomer. |

Use `charger: nagl`. {term}`NAGL` predicts {term}`AM1-BCC` charges and needs
no extra program. The `am1bcc` method needs AmberTools, and `espaloma` needs
espaloma-charge. The PolyzyMD environments include neither.

For every other key of the config, see {doc}`../reference/configuration`.

## Build and run the system

1. Validate the config:

   ```bash
   polyzymd validate -c config.yaml
   ```

2. Build one replicate. The first build generates the fragments and the chains,
   and can take several minutes:

   ```bash
   polyzymd build -c config.yaml -r 1
   ```

3. Run the simulation. On a SLURM cluster, submit it from the `build`
   environment:

   ```bash
   polyzymd submit -c config.yaml -r 1 --preset aa100
   ```

   The jobs start the simulation environment themselves. See
   {doc}`hpc_slurm`. To run on your own computer, use
   `polyzymd run -c config.yaml -r 1` instead.

## What the build does

For the monomers, the build does these steps once:

1. It reads the monomer SMILES.
2. It applies the initiation template to each monomer. For ATRP, this adds a
   chlorine at the double bond.
3. It applies the polymerization template to make the middle fragments, which
   have two bonding sites.
4. It applies the termination template to make the end fragments, which have
   one bonding site and keep the double bond.
5. It saves the fragments as `<type_prefix>_monomer_group.json` in
   `cache_directory`.

For each chain, the build does these steps:

1. It draws a sequence from the monomer probabilities. The replicate number
   seeds the draw, unless you set `random_seed`.
2. It puts end fragments at the two ends and middle fragments between them.
3. It builds the 3D structure with Polymerist. The sequence seeds the
   structure, so the same sequence always gives the same coordinates.
4. It checks that no bond passes through a ring. If one does, it builds the
   chain again, up to `max_retries` times.
5. It assigns partial charges with `charger`.
6. It saves the charged chain as
   `<type_prefix>_seq=<sequence>_<length>-mer_charged.sdf` in
   `cache_directory`, with a `.metadata.json` file beside it. This is the name
   that `cached` mode reads, so a later config can set
   `generation_mode: cached` and `sdf_directory` to this folder.

The next build with the same monomers loads the fragments and chains from
`cache_directory`. It generates a chain again if the metadata does not match
the config.

## Use another monomer

The bundled templates are for methacrylates. To add a methacrylate monomer:

1. Write its SMILES with the C=C double bond, and with explicit hydrogens.
2. Add it to `monomers` with a new `label`:

   ```yaml
   monomers:
     - label: "C"
       probability: 0.1
       name: "MyMonomer"
       smiles: "[H]C([H])=C(C(=O)O...)C([H])([H])[H]"
   ```

3. Adjust the other probabilities so that all of them sum to 1.0.

For another chemistry, such as ring-opening polymerization, write your own
three `.rxn` templates and give their paths in `reactions`.

## Fix common problems

### "No 1-site terminal fragment found for monomer"

The fragment cache does not match the monomers of the config. Delete the cache
and build again:

```bash
rm -rf .polymer_cache
```

### "No monomer name configured for sequence label"

A monomer has no `name`. Give every monomer a `name`. `polyzymd validate`
passes such a config, so this error first appears at build time.

### "Failed to build polymer after N attempts due to ring-piercing"

1. Raise `max_retries`.
2. Use shorter chains.
3. Check that each monomer SMILES is correct.

### "Dynamic generation mode requires 'length' >= 3"

Set `length` to 3 or more.

## See also

- {doc}`polymers`: the keys that both modes share, and how the chains are
  placed.
- {doc}`../reference/configuration`: every key of the config.
- {doc}`run_gromacs`: run the system with GROMACS.
