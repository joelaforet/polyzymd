# Add polymers to a simulation

Add random co-polymer chains around the protein with a `type: polymer` entry
in the `solutes:` list of `config.yaml`. The builder draws a monomer sequence for each chain, places
the chains around the protein with {term}`PACKMOL`, and then adds water and
ions. No SDF files of your chains? Use `generation_mode: dynamic`, which
builds each chain from monomer SMILES: see {doc}`dynamic_polymers`. This page
covers the `cached` mode, which reads pre-built SDF files.

:::{admonition} Environment setup
:class: tip

Run the commands on this page in the `build` environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

PolyzyMD gets the polymer structures in one of two ways. The
`generation_mode` key selects the way:

| `generation_mode` | Source of each chain | Guide |
|---|---|---|
| `cached` (default) | A pre-built SDF file for each sequence, in `sdf_directory` | This page |
| `dynamic` | Built with mBuild from star-sited monomer and cap SMILES | {doc}`dynamic_polymers` |

## Write the polymer entry

This example adds two 5-mer chains of SBMA and EGPMA from pre-built SDF files:

```yaml
solutes:
  # ... the protein entry ...
  - type: polymer
    type_prefix: "SBMA-EGPMA"
    monomers:
      - label: "A"
        probability: 0.98
        name: "SBMA"    # sulfobetaine methacrylate
      - label: "B"
        probability: 0.02
        name: "EGPMA"   # ethylene glycol phenyl ether methacrylate
    length: 5           # monomers per chain
    count: 2            # number of chains
    sdf_directory: "polymer_sdfs/SBMA-EGPMA"
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: preset}   # the charges stored in the SDF files
```

The entry is one item of `solutes:`; all its keys are indented under the `-`.

The keys have these meanings:

| Key | Meaning |
|---|---|
| `type_prefix` | The polymer name. It is the first part of each SDF file name. |
| `monomers` | One entry for each monomer type. `label` is one character. |
| `probability` | The chance that a position in the chain gets this monomer. The probabilities must sum to 1.0. |
| `name` | The monomer name. Its first three letters, in upper case, become the residue name unless you set `residue_name`. |
| `length` | The number of monomers in each chain. |
| `count` | The number of chains in the box. |
| `sdf_directory` | The folder with the pre-built SDF files. `cached` mode requires it. |
| `force_field` | The force field file, as a list. It is the same file as `solvent.force_field`. |
| `charges` | `{method: preset}` in `cached` mode: the builder uses the charges in the SDF files. |

A relative `sdf_directory` is relative to the folder that holds `config.yaml`.

For a homopolymer, give one monomer with `probability: 1.0`.

## Name the SDF files

For each chain, the builder draws a sequence such as `AABAA`. In `cached`
mode it looks up a sequence and its reverse under one file name: the form that
comes first in alphabetical order, so `ABAAA` is looked up as `AAABA`. If that
file is missing, it loads the file of the reversed sequence, which is how
`dynamic` mode names a chain it built in that direction.

For a head-to-tail chain a sequence and its reverse are different molecules
(the head sits on a different monomer), so in `cached` mode the file you
provide sets the direction of every chain of that sequence. The two are
equivalent only in distribution, and only when the head and tail ends are
drawn alike. `dynamic` mode reads each sequence from head to tail and builds
it as drawn; see {doc}`dynamic_polymers`.

When the config gives `end_groups`, `cached` mode also draws each chain's
caps, and the sequence starts with the head cap's label and ends with the tail
cap's label, as in `dynamic` mode (`length` then counts the two ends).

The builder then loads this file from `sdf_directory`:

```text
<type_prefix>_seq=<sequence>_<length>-mer_charged.sdf
```

For the example above, the file of sequence `AAAAA` is
`SBMA-EGPMA_seq=AAAAA_5-mer_charged.sdf`. If the file is not in
`sdf_directory`, the builder looks in `cache_directory` (default
`.polymer_cache`). If the file is in neither folder, the build stops with
`FileNotFoundError` and names both paths. `cached` mode does not generate a
missing chain.

Provide a file for every sequence that the probabilities can produce. The
builder uses the partial charges stored in each SDF file.

## Control the sequence draw

The replicate number seeds the sequence draw, so each replicate gets its own
set of chains. To give every replicate the same chains, set a fixed seed:

```yaml
solutes:
  - type: polymer
    # ... the other polymer keys ...
    random_seed: 42
```

The replicate number still seeds the PACKMOL placement.

## Build a system without polymers

For a control without polymer, leave the polymer entry out of `solutes:`.

## How the builder places the chains

The builder does these steps:

1. It computes the periodic cell from the protein and substrate. The edge is
   their diameter plus 2 x (`packing.padding` of the polymer solute +
   `solvent.box.padding`). Expect a much larger box than without polymers:
   with the default 2 nm polymer padding, the quickstart's Trp-cage in water
   (6,138 atoms) became 42,323 atoms with two SBMA-EGMA 5-mers.
2. It centers the bounding box of the protein and substrate in the rectangular
   brick of that cell. It holds them fixed.
3. It packs the chains inside the brick and inside a sphere around the solute.
   PACKMOL keeps each polymer atom at least `packing.tolerance`
   (default 2.0 Å) from the solute and from other chains.
4. It fills the remaining space of the brick with water. It does not move the
   system between packing and solvation.
5. It adds ions to neutralize the system and to reach the salt concentration.

The replicate number seeds PACKMOL, so each replicate gets its own placement.
The cell does not depend on the placement. All replicates of one condition
therefore get the same box, the same number of waters and the same number of
ions.

### Packing settings

These keys go under `packing:` in the polymer entry:

```yaml
solutes:
  - type: polymer
    # ... the other polymer keys ...
    packing:
      padding: 2.0            # nm of room for the chains (default 2.0)
      tolerance: 2.0          # Å, smallest distance between molecules (default 2.0)
      confine_to_sphere: true # keep the chains in a sphere around the solute (default)
      movebadrandom: false    # PACKMOL movebadrandom keyword (default false)
      nloop: 200              # PACKMOL loops per molecule type (default 200)
```

The sphere has a radius of half the solute bounding-box diagonal plus
`padding`. With `confine_to_sphere: false`, the chains fill the whole brick.
In the SBMA-EGMA pentamer systems, the mean polymer-to-protein distance was
then 23 Å instead of 19 Å.

Packing inside the brick keeps each chain away from its own periodic images.
Two atoms inside the brick, shrunk by `tolerance`, are at least `tolerance`
apart across every lattice vector.

To give the chains more room, raise `padding`. The cell is always computed
before packing, so every replicate of a condition gets the same cell.

### Checks after packing

PACKMOL can exit with code 173 ("imperfect packing") for a dense system. The
build accepts the result with a warning, and energy minimization removes the
remaining close contacts. These checks run before any simulation:

- The build stops with `SolvationClashError` if more than 20 packed atoms lie
  within half the tolerance of the solute. This many contacts means that the
  solute and the packed molecules are in different coordinate frames.
- The build stops with `PeriodicImageClashError` if an atom lies within half
  the tolerance of one of its own periodic images. The build checks this after
  packing and again after solvation, over all 26 neighbor cells.
- Minimization holds the protein and substrate heavy atoms fixed by default.
  Water and polymer atoms move to remove the close contacts.

See {doc}`../explanation/simulation_safeguards` for the reasons behind these
checks.

## Check the system before you run it

1. Validate the config:

   ```bash
   polyzymd validate -c config.yaml
   ```

2. Preview the build. The preview lists the polymer and what the replicate
   number seeds:

   ```bash
   polyzymd build -c config.yaml --dry-run
   ```

3. Build one replicate:

   ```bash
   polyzymd build -c config.yaml -r 1
   ```

4. Open `solvated_system.pdb` in the {term}`replicate folder`. Polymer atoms
   are chain `C`.

A system with more polymer atoms runs more slowly. Test a small system, such as
two 5-mers, before you build a large one.

## Fix common problems

### PACKMOL does not find a solution

Give the chains more room, or let PACKMOL work longer:

```yaml
solutes:
  - type: polymer
    # ... the other polymer keys ...
    packing:
      padding: 2.5          # more room for the chains
      movebadrandom: true   # helps when there are many different chain sequences
      nloop: 400            # raise when PACKMOL stops at the loop limit
solvent:
  box:
    padding: 2.0          # a larger box
```

You can also lower `count` or `length`.

An exit code of 173 is not a failure. The build continues. A
`SolvationClashError` is a different problem; see {doc}`troubleshooting`.

### OpenFF cannot assign parameters to a chain

1. Check that each SDF file has correct bond orders, formal charges and
   explicit hydrogens.
2. Check that each SDF file has partial charges. The builder uses them.
3. Build the chains in `dynamic` mode with
   `charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}`, so that
   {term}`NAGL` assigns the charges. NAGL needs no extra program. See
   {doc}`dynamic_polymers`.

The `am1bcc` charge method ({term}`AM1-BCC`) needs AmberTools, which the
PolyzyMD environments do not include. Use it only in an environment where you
installed AmberTools.

### The simulation is unstable with polymers

1. Add a stage that holds the protein while the polymer relaxes. See
   {doc}`equilibration`.
2. Make the free equilibration stage longer.
3. Open `solvated_system.pdb` and look for chains that pass through the protein
   or through a ring.

## See also

- {doc}`dynamic_polymers`: build the chains from monomer SMILES.
- {doc}`equilibration`: hold the protein while the polymer relaxes.
- {doc}`../reference/configuration`: every key of the polymer solute.
