# Author Branched Polymers Offline with the CGSmiles Notebook

Use this guide when you need a branched, tree-like, or dendrimer-style polymer
that PolyzyMD should pack as a pre-authored molecule. The notebook created by
`polyzymd init` is an offline authoring helper: it writes a charged SDF and a
`provided_molecules` YAML snippet. Runtime simulation builds do not execute
CGSmiles.

## Create the scaffold

Start from a normal PolyzyMD project:

```bash
pixi run -e build polyzymd init --name branched_polymer_project
cd branched_polymer_project
```

The scaffold includes:

```text
branched_polymer_project/
├── config.yaml
├── structures/
├── notebooks/
│   └── cgsmiles_polymer_scaffold.ipynb
├── generated_molecules/
├── job_scripts/
└── slurm_logs/
```

The scaffolded `config.yaml` is intentionally incomplete. Fill in the required
enzyme, thermodynamics, solvent, simulation, and polymer sections before running
validation or build commands.

Open `notebooks/cgsmiles_polymer_scaffold.ipynb` in the build environment or a
compatible Jupyter environment. The current build environment includes the
notebook stack, CGSmiles 1.0.2, mBuild development pin, HOOMD, and OpenFF tools
used by the scaffold.

Exact launch command from the repository or any shell with pixi available:

```bash
pixi run -e build jupyter lab branched_polymer_project/notebooks/
```

## Edit only the three user cells

The notebook has five cells. Edit the three user-facing cells:

1. **Graph cell** - define the offline mBuild coarse graph, such as a tree,
   star, or dendrimer-like graph.
2. **Fragments/templates/maps cell** - provide CGSmiles fragments, atomistic
   templates, exact atom maps, port maps, and leaving maps.
3. **Output/build cell** - choose the output name and run the single build call.

The notebook does not suggest chemistry. You must provide exact atom, port, and
leaving maps for the polymer you are authoring.

The output cell resolves paths from the project root and writes under
`generated_molecules/`. Output writes are atomic so an interrupted run does not
leave a partially written SDF in place. The notebook stays deliberately simple;
if you need a specialized RDKit or OpenMM relaxation workflow, do that as a
separate validation step before using the SDF in production.

## Stay within current support

Current notebook support is intentionally narrow:

- each coarse edge resolves to exactly one atomistic interfragment record;
- edge `order` is bond order and may be `1`, `1.5` aromatic, `2`, or `3`;
- multiple distinct interfragment records, PIM-style multi-point residue pairs,
  and ring closures/cyclic graphs are not supported by the scaffold.

For unsupported cases, author and charge the molecule with an external workflow,
then provide the charged single-molecule SDF through `provided_molecules`.

## Use the generated molecule in `config.yaml`

After the notebook succeeds, copy or paste the emitted YAML into your polymer
configuration. The simulation runtime consumes the charged SDF; it does not
re-run CGSmiles.

For provided-only systems, paste the full emitted block:

```yaml
polymers:
  enabled: true
  generation_mode: "provided"
  provided_molecules:
    - name: "authored_branched_polymer"
      entries:
        - sdf_path: "generated_molecules/authored_branched_polymer.sdf"
          count: 1
```

If you already use `generation_mode: "dynamic"`, `"fragments"`, or `"cached"`,
keep the generated-polymer fields in that existing block and merge only the
`provided_molecules` list from the notebook output.

Then validate and build as usual:

```bash
pixi run -e build polyzymd validate -c config.yaml
pixi run -e build polyzymd build -c config.yaml -r 1
```

## Validation boundary

The authored SDF is validated through OpenFF before use. It must be one
connected molecule with finite coordinates, no dummy atoms, and complete finite
partial charges. See {doc}`../reference/data_requirements` for the full SDF
requirements and {doc}`../reference/polymer_mbuild_openff` for the adapter
boundary.

Current end-to-end validation covered one authored molecule through OpenFF
parameterization, Packmol packing/solvation with TIP3P, and completed a
1000-step CPU OpenMM smoke run. Treat that as evidence for the documented path,
not as general support for arbitrary cyclic, PIM, or multi-record chemistry.

## See also

- {doc}`polymers` - add generated and provided polymers to simulations
- {doc}`dynamic_polymers` - generate native linear methacrylate polymers
- {doc}`../reference/configuration` - polymer configuration fields
- {doc}`../reference/polymer_mbuild_openff` - mBuild/OpenFF adapter reference
