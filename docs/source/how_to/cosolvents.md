# Add a co-solvent

Add co-solvent molecules, such as urea or DMSO, to the water of a condition.
The usual way is to copy the water condition of a study and add a
`co_solvents:` entry to the `solvent:` section of the copy. PACKMOL then
places the co-solvent with the water and the ions. This works on OpenMM and
GROMACS.

```{figure} ../_static/renders/cosolvent.webp
:figclass: pz-render
:class: dark-light
:alt: Trp-cage drawn as a gray cartoon in a box of water, shown as small light spheres, with urea molecules with green carbons and sodium and chloride ions as purple and yellow spheres.

The shipped example `examples/cosolvent`: Trp-cage (gray) in 2 M urea (green
carbons) and water, with Na⁺ (purple) and Cl⁻ (yellow). Water is drawn as its
oxygens; the water and urea in front of the protein are cut away.
```

:::{admonition} Environment setup
:class: tip

Run the commands on this page in the `build` environment:

```bash
pixi shell -e build
```
:::

## Copy the condition

Make the co-solvent condition a copy of the water condition, so that the two
differ only in the co-solvent:

```bash
polyzymd study add-condition "Urea 2 M" --from Water --study trpcage
```

The copy is `trpcage/conditions/urea_2_m/config.yaml`. Change its
`description:`.

## Write the entry

Add the entry to the `solvent:` section, at the same indent as `water:`. If
the section has a line `co_solvents: []`, replace that line:

```yaml
solvent:
  force_field: ["openff-2.3.0.offxml"]
  water:
    model: "tip3p"
  co_solvents:
    - name: "urea"
      concentration: 2.0          # mol/L
      charges: {method: preset}   # the charges PolyzyMD ships for urea
  ions:
    neutralize: true
    nacl_concentration: 0.15
```

Each co-solvent needs three things:

| Key | What to write |
|---|---|
| `name` | A name from the built-in library, such as `urea`, `dmso`, `ethanol` or `glycerol`, or your own name together with `smiles:` |
| one amount | `mole_fraction` (0 to 1, replaces water), `concentration` (mol/L, added on top of the water) or `count` (molecules, added on top of the water). Give exactly one |
| `charges` | `{method: preset}` for a library co-solvent. For your own molecule, `{method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}`, `{method: am1bcc}` or `{method: espaloma}` |

For example, 10 mol% DMSO and a molecule PolyzyMD does not ship:

```yaml
  co_solvents:
    - name: "dmso"
      mole_fraction: 0.10
      charges: {method: preset}
    - name: "ethyl_acetate"
      smiles: "CCOC(=O)C"
      mole_fraction: 0.05
      charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

The co-solvent is parameterized with the force field file of
`solvent.force_field`. For the library, the formulas behind each amount and
the residue names, see [Co-solvents](../reference/configuration.md#co-solvents) in the
configuration reference.

## Check and build

```bash
polyzymd validate -c trpcage/conditions/urea_2_m/config.yaml
polyzymd build -c trpcage/conditions/urea_2_m/config.yaml -r 1 --dry-run
```

The summary of `validate` gives `Co-solvents: urea`. The dry run prints the
residue name of each co-solvent; urea is `URE`. Then run the replicates as
for any condition.

The build log says how many molecules it added. For the quickstart box of
Trp-cage, 2 M urea is:

```
Adding 1941 water, 5 Na+, 6 Cl-
Adding 75 urea molecules (2.0 M)
```

`build_provenance.md` lists the co-solvent among the components. Its physics
audit reports the density as `INFO`. With `concentration` or `count`, the
density is above the target of `box.target_density`, because the co-solvent
is added on top of the water; about 1.12 g/mL for 2 M urea. This is
expected.

## Analyze

Select the co-solvent by its residue name, such as `resname URE`. For
example, use it in a SASA context ({doc}`analysis_sasa_quickstart`) or as a
group of `hydrogen_bonds` ({doc}`hydrogen_bonds`). The study compares the
co-solvent condition with the water control, as for a polymer condition.

## See also

- {doc}`../tutorials/own_system`: a study with a water and a 2 M urea
  condition, from setup to analysis.
- {ref}`the shipped examples <shipped-examples>`: the shipped co-solvent example.
- [Co-solvents](../reference/configuration.md#co-solvents) in the configuration reference: every key,
  the library and the formulas.
