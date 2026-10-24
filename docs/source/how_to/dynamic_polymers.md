# Generate polymers from SMILES

Build polymer chains from monomer SMILES with `generation_mode: "dynamic"`.
In this mode you do not need an SDF file for each sequence. PolyzyMD builds
each chain with mBuild from the repeat unit of each monomer and the two caps,
and assigns partial charges with the method in `charges`.

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

## Write the monomers as star-sited SMILES

Give each monomer as its **repeat unit**, as it is in the chain, with `*` on
the two atoms that bond to the neighbouring repeat units. For a vinyl,
acrylate or methacrylate monomer these are the two carbons of the C=C double
bond, which is a single bond in the chain:

| Monomer | SMILES as sold (PubChem) | Repeat unit |
|---|---|---|
| methyl methacrylate (MMA) | `C=C(C)C(=O)OC` | `[*]CC([*])(C)C(=O)OC` |
| 2-hydroxyethyl methacrylate (HEMA) | `C=C(C)C(=O)OCCO` | `[*]CC([*])(C)C(=O)OCCO` |
| sulfobetaine methacrylate (SBMA) | `C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]` | `[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]` |

The first `*` bonds toward the head of the chain, the second toward the tail.
Charges in the SMILES are kept: an SBMA repeat unit keeps its N+ and O-.

The notebook {doc}`../tutorials/star_sited_monomers` draws a monomer with its
atom indices and puts the `*` on the atoms you name.

## Choose the caps

Each chain has two ends. `end_groups` gives the cap of each end as SMILES with
one `*` on the atom that bonds to the chain. There is no default: the caps are
part of the molecule you simulate.

| Cap | SMILES |
|---|---|
| hydrogen | `[*][H]` |
| methyl | `[*]C` |

`head` caps the first repeat unit of a sequence, `tail` the last.

## Write the polymer entry

```yaml
solutes:
  # ... the protein entry ...
  - type: polymer
    generation_mode: "dynamic"
    type_prefix: "SBMA-EGPMA"
    monomers:
      - label: "A"
        probability: 0.7
        name: "SBMA"
        residue_name: "SBM"
        smiles: "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
      - label: "B"
        probability: 0.3
        name: "EGPMA"
        residue_name: "EGP"
        smiles: "[*]CC([*])(C)C(=O)OCCOc1ccccc1"
    end_groups:
      head: "[*][H]"
      tail: "[*][H]"
    length: 5
    count: 2
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
    max_retries: 10
    cache_directory: ".polymer_cache"
```

`dynamic` mode adds these keys and rules:

| Key | Meaning |
|---|---|
| `smiles` | The repeat unit of the monomer, with exactly two `*`. Every monomer needs one. |
| `end_groups` | `head` and `tail`, each a cap with exactly one `*`. Required. |
| `name` | The monomer name, used in log messages and the cache metadata. Optional; the `label` is used without it. |
| `residue_name` | The 3-character residue name of the monomer in the topology. Optional. If you leave it out, PolyzyMD makes one from `name`. It must differ between monomers and must not be an amino-acid, nucleic-acid, water or ion name. |
| `charges` | The partial-charge method for the chains: `{method: nagl, model: ...}`, `{method: am1bcc}` or `{method: espaloma}`. `preset` is for `cached` mode only. |
| `max_retries` | The number of builds tried, each with new monomer conformers, before a chain with a pierced ring or overlapping atoms is refused (default 10). |
| `cache_directory` | The folder of the generated chains (default `.polymer_cache`, relative to the config). |

{term}`NAGL` predicts {term}`AM1-BCC` charges and needs no extra program.
Its `model` is required, such as `openff-gnn-am1bcc-1.0.0.pt`, the model
Sage 2.3.0 names. The `am1bcc` method needs AmberTools, and `espaloma` needs
espaloma-charge. The PolyzyMD environments include neither.

For every other key of the config, see {doc}`../reference/configuration`.

## Build and run the system

1. Validate the config. This reads every monomer and cap SMILES with RDKit,
   checks the number of `*`, and checks that the charge method covers every
   element:

   ```bash
   polyzymd validate -c config.yaml
   ```

2. Build one replicate. Each new chain takes a few seconds:

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

For each chain, the build does these steps:

1. It draws a sequence from the monomer probabilities. The replicate number
   seeds the draw, unless you set `random_seed`. The sequence is read from
   head to tail, so `AAB` and `BAA` are different chains.
2. It joins the repeat units of the sequence and the two caps: each bond
   between two units replaces the `*` of both. The bond orders and formal
   charges are those of your SMILES.
3. It builds the 3D structure with mBuild's `Polymer`. The sequence seeds the
   structure, so the same sequence always gives the same coordinates.
4. It assigns partial charges with the method in `charges`, then relaxes the
   chain in vacuum with its `force_field` (an energy minimization).
5. It checks that no bond passes through a ring and that no two atoms more
   than two bonds apart are closer than 1.0 Å. If a check fails, it builds the
   chain again with new monomer conformers, up to `max_retries` times.
6. It saves the charged chain as
   `<type_prefix>_seq=<sequence>_<length>-mer_charged.sdf` in
   `cache_directory`, with a `.metadata.json` file beside it. This is the name
   that `cached` mode reads, so a later config can set
   `generation_mode: cached`, `sdf_directory` to this folder and
   `charges: {method: preset}`.

Each repeat unit is one residue named by its `residue_name`. The head cap
belongs to the first residue and the tail cap to the last. Repeat units join
head to tail: the second `*` of one unit bonds to the first `*` of the next.

Tacticity is not controlled. Every repeat unit of one monomer is built from
the same 3D structure, so its backbone stereocentres all get the same
configuration: a chain of one methacrylate comes out isotactic, with the
handedness set by the sequence's seed. The chain SDF records the
stereochemistry of every atom.

The next build with the same monomers loads the chains from
`cache_directory`. It builds a chain again if the metadata does not match the
config, for example after you change a SMILES, a cap, the charges or the
mBuild version.

The build record (`build_provenance.json` and `.md`) names the mBuild version
and commit that built the chains, and lists mBuild under "Methods to cite".

## Fix common problems

### "has 0 '*' sites; a repeat unit needs exactly 2"

The monomer SMILES is the monomer as sold, not the repeat unit. Put `*` on the
two linking atoms; see {doc}`../tutorials/star_sited_monomers`.

### "both '*' ... are on one atom"

Each `*` must be on a different atom. Write a longer repeat unit, such as
`[*]CC[*]` for polyethylene.

### "no build out of N passed the geometry checks"

1. Raise `max_retries`.
2. Check each monomer SMILES in the notebook.

### "dynamic polymers need mBuild from https://github.com/joelaforet/mbuild"

The environment has another mBuild. Run PolyzyMD in its `build` pixi
environment, which installs the right one.

## See also

- {doc}`polymers`: the keys that both modes share, and how the chains are
  placed.
- {doc}`../tutorials/star_sited_monomers`: write star-sited SMILES with a
  notebook.
- {doc}`../reference/configuration`: every key of the config.
- {doc}`run_gromacs`: run the system with GROMACS.
