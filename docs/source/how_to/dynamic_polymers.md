# Generate polymers from SMILES

Build polymer chains from monomer SMILES with `generation_mode: "dynamic"`.
In this mode you do not need an SDF file for each sequence. PolyzyMD builds
each chain with mBuild from the repeat unit of each monomer and the two caps,
and assigns partial charges with the method in `charges`.

```{figure} ../_static/renders/polymer_system.webp
:figclass: pz-render
:class: dark-light
:alt: Two coiled polymer chains drawn as thick blue sticks with red, blue and yellow sulfonate and ammonium groups, beside a small protein drawn as a gray cartoon, in water shown as small light spheres with a few purple sodium and yellow chloride ions.

The shipped example `examples/polymer`: two SBMA 5-mers built from star-sited
SMILES (carbons blue; the sulfonate S yellow and O red, the ammonium N blue)
beside Trp-cage (gray), in water with Na⁺ (purple) and Cl⁻ (yellow). These are
the coordinates as built, before any simulation; water and ions in front of
the chains and the protein are cut away.
```

For the keys that both modes share, such as `monomers`, `length`, `count` and
`packing`, see {doc}`polymers`.

:::{admonition} Environment setup
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

Each chain has two ends. `end_groups` lists the caps each end can get, as
SMILES with one `*` on the atom that bonds to the chain, and the probability
of each. Every chain draws its head cap and its tail cap from these lists. The
probabilities of each end sum to 1. A fixed cap is one option with probability
1.0. There is no default: the caps are part of the molecule you simulate.

| Cap | Option |
|---|---|
| hydrogen | `{label: H, smiles: "[*][H]", probability: 1.0}` |
| methyl | `{label: M, smiles: "[*]C", probability: 1.0}` |
| an SBMA unit that kept its C=C | `{label: A, smiles: "[*]C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]", residue_name: SBM, probability: 0.7}` |

`label` is the one character the chain's sequence shows for that cap. A cap
with a `residue_name` is its own residue; a cap without one joins the residue
of the repeat unit it caps.

`length` counts the two ends: `length: 5` is a head cap, three repeat units and
a tail cap. This holds for every cap, a hydrogen included, so `length` is at
least 3.

### Chain ends that keep the monomer's double bond

Before 1.3.0, PolyzyMD built dynamic chains with ATRP reaction templates.
Those chains, such as the SBMA-EGMA 5-mers of studies made with those
releases, are three polymerized repeat units with a whole monomer at each end
that kept its C=C double bond, bonded to the chain through its CH2 carbon (now
`=CH-`), and each end's monomer was drawn with the monomer probabilities. This
is a modelling choice: it stands for a termination step that restores the
double bond. A real ATRP chain ends in an initiator fragment and a halide. In
star-sited SMILES the end form of a methacrylate is `[*]C=C(C)C(=O)OR`, and
the entry below builds those 5-mers. EGMA is methoxy tetra(ethylene glycol)
methacrylate.

## Write the polymer entry

```yaml
solutes:
  # ... the protein entry ...
  - type: polymer
    generation_mode: "dynamic"
    type_prefix: "SBMA-EGMA"
    monomers:
      - label: "A"
        probability: 0.7
        name: "SBMA"
        residue_name: "SBM"
        smiles: "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
      - label: "B"
        probability: 0.3
        name: "EGMA"
        residue_name: "EGM"
        smiles: "[*]CC([*])(C)C(=O)OCCOCCOCCOCCOC"
    end_groups:
      head:
        - {label: "A", name: "SBMA", residue_name: "SBM", smiles: "[*]C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]", probability: 0.7}
        - {label: "B", name: "EGMA", residue_name: "EGM", smiles: "[*]C=C(C)C(=O)OCCOCCOCCOCCOC", probability: 0.3}
      tail:
        - {label: "A", name: "SBMA", residue_name: "SBM", smiles: "[*]C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]", probability: 0.7}
        - {label: "B", name: "EGMA", residue_name: "EGM", smiles: "[*]C=C(C)C(=O)OCCOCCOCCOCCOC", probability: 0.3}
    length: 5
    count: 2
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
    max_retries: 10
    cache_directory: ".polymer_cache"
```

A chain drawn as `AABAB` has an SBMA head end, the repeat units SBMA, EGMA,
SBMA, and an EGMA tail end: five residues, SBM SBM EGM SBM EGM. Its file is
`SBMA-EGMA_seq=AABAB_5-mer_charged.sdf`.

For hydrogen caps on both ends, write:

```yaml
    end_groups:
      head: [{label: "H", smiles: "[*][H]", probability: 1.0}]
      tail: [{label: "H", smiles: "[*][H]", probability: 1.0}]
```

`dynamic` mode adds these keys and rules:

| Key | Meaning |
|---|---|
| `smiles` | The repeat unit of the monomer, with exactly two `*`. Every monomer needs one. |
| `end_groups` | `head` and `tail`, each a list of caps with `label`, `smiles` (exactly one `*`), `probability`, and optionally `name` and `residue_name`. Each end's probabilities sum to 1. Required. |
| `length` | Monomers per chain, the two ends included; at least 3. |
| `name` | The monomer name, used in log messages and the cache metadata. Optional; the `label` is used without it. Names must differ between monomers. |
| `residue_name` | The 3-character residue name of the monomer in the topology. Optional. If you leave it out, PolyzyMD makes one from `name`. It must not be an amino-acid, nucleic-acid, water or ion name, and it must differ between monomers, except between stereoisomers of one monomer (see {ref}`tacticity`). |
| `charges` | The partial-charge method for the chains: `{method: nagl, model: ...}`, `{method: am1bcc}` or `{method: espaloma}`. `preset` is for `cached` mode only. |
| `max_retries` | The number of builds tried, each with new monomer conformers, before a chain with a pierced ring or overlapping atoms is refused (default 10). |
| `cache_directory` | The folder of the generated chains (default `.polymer_cache`, relative to the config). |

{term}`NAGL` predicts {term}`AM1-BCC` charges and needs no extra program.
Its `model` is required, such as `openff-gnn-am1bcc-1.0.0.pt`, the model
Sage 2.3.0 names. The `am1bcc` method needs AmberTools, and `espaloma` needs
espaloma-charge. The PolyzyMD environments include neither.

For every other key of the config, see {doc}`../reference/configuration`.

(tacticity)=
## Set the tacticity

A stereocentre written in the repeat unit's SMILES is kept in the chain. `@`
and `@@` fix the configuration of the backbone carbon relative to its
neighbours as written, so two units with the same SMILES have the same
relative configuration (a meso dyad) and two with opposite marks the opposite
(a racemo dyad). The CIP label (R or S) of a carbon in the chain can differ
from that of the unit, because it depends on the chain around it; tacticity is
about relative configuration only. To control it, give each configuration of
the backbone stereocentre its own label, with the same `residue_name`:

```yaml
    monomers:
      - {label: "X", probability: 0.5, name: "MMA-x", residue_name: "MMA", smiles: "[*]C[C@@]([*])(C)C(=O)OC"}
      - {label: "Y", probability: 0.5, name: "MMA-y", residue_name: "MMA", smiles: "[*]C[C@]([*])(C)C(=O)OC"}
```

Each repeat unit then draws its configuration independently: equal
probabilities give an atactic chain, and one label with probability 1 an
isotactic chain. A syndiotactic chain (strictly alternating) cannot be drawn
this way; build those chains yourself and give them in `cached` mode.

:::{warning}
A repeat unit written without its stereocentre gets one configuration from the
builder: within one chain every unit of that monomer has the same
configuration, but two monomers, or the same monomer in two chains, can differ.
The tacticity of such chains is not controlled; give stereo letters to control
it.
:::

The same holds for a double bond: written with `/` and `\` (for example
`[*]/C=C(/C)C(=O)OC`, the Z end), its configuration is kept; written without
them (`[*]C=C(C)C(=O)OC`), the builder picks one. Each chain's SDF records the
configuration of every stereocentre and double bond as built.

Each chain is drawn monomer by monomer, so the build makes random (statistical)
copolymers. To simulate block copolymers, build the chains yourself and give
them in `cached` mode, or open an issue on GitHub so we can add them.

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

1. It draws the chain's head cap, `length - 2` repeat units and its tail cap
   from their probabilities. The replicate number seeds the draw, unless you
   set `random_seed`. The sequence is the head cap's label, the repeat units'
   labels and the tail cap's label, read from head to tail, so `AAB...` and
   `BAA...` are different chains.
2. It joins the units: each bond between two units replaces the `*` of both.
   Repeat units join head to tail: the second `*` of one unit bonds to the
   first `*` of the next. The bond orders and formal charges are those of your
   SMILES.
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
   `charges: {method: preset}`; cached mode then checks the config's monomer
   SMILES and caps, if it gives them, against that metadata.

Each repeat unit is one residue named by its `residue_name`; a cap is its own
residue when it has a `residue_name`. The chain SDF records the configuration
of every stereocentre and double bond as built; a build whose relaxation would
change one is refused and built again.

The next build with the same monomers loads the chains from
`cache_directory`. It builds a chain again if the metadata does not match the
config, for example after you change a SMILES, a cap, the charges or the
mBuild version.

The build record (`build_provenance.json` and `.md`) names the mBuild version
and commit that built the chains, and lists mBuild under "Methods to cite".

## Fix common problems

### "end_groups.head probabilities must sum to 1.0"

Each end's caps are drawn from their probabilities. Make them sum to 1; a
single cap needs `probability: 1.0`.

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
