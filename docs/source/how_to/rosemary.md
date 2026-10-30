# Parameterize a protein with OpenFF Rosemary alpha

```{index} !rosemary, !openff-3.0.0, !protein force field
```

Parameterize a whole protein, with its mutations and modifications, with
OpenFF's protein force field in development instead of ff14SB, in OPC3 water.
The fields are listed in [OpenFF Rosemary alpha](protein-rosemary-alpha) of
the configuration reference.

```{warning}
Rosemary is an alpha release, not for production. The audit reports it
(`parameters.rosemary_alpha`).
```

:::{note}
Run `polyzymd` from the `build` pixi environment (`pixi run -e build ...`).
:::

## Step 1: Use Rosemary for the protein

```yaml
solutes:
  - type: protein
    pdb_path: "structures/enzyme.pdb"
    force_field: ["openff_no_water-3.0.0-alpha2b.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
```

`alpha2b` is the newest candidate the pixi environments install; earlier
candidates such as `alpha0` are still accepted. Rosemary must be the
protein's only force field file. NAGL charges the whole protein molecule.

## Step 2: Use OPC3 water

OpenFF pairs Rosemary with OPC3 water. Set the water model:

```yaml
solvent:
  force_field: ["openff-2.3.0.offxml"]
  water:
    model: "opc3"
```

PolyzyMD loads the force fields in this order: the protein's (Rosemary), then
`solvent.force_field`, then `opc3.offxml`. For the protein and the water this
is `ForceField("openff_no_water-3.0.0-alpha2b.offxml", "opc3.offxml")`. The
water and the Na+ and Cl- ions take their parameters from `opc3.offxml` (OPC3
water and the Li/Merz 12-6 ions fitted for it); see
[Water Models](water-models).

With any other water model the build stops at the audit:

```text
FAIL  solvent.rosemary_water  system  trpcage (openff_no_water-3.0.0-alpha2b.offxml) is on
      OpenFF Rosemary, which OpenFF pairs with OPC3 water, but the water is tip3p. ...
```

To keep another water model on purpose, for example to compare with earlier
TIP3P runs, accept the check with a reason:

```yaml
modeling_exceptions:
  - check: solvent.rosemary_water
    reason: "TIP3P on purpose: compared with our earlier TIP3P runs of this enzyme"
```

The build record then lists the finding as ACCEPTED, with the reason.

## Step 3: Add changes without parameter sets

Mutations and modifications work as on ff14SB, but the protein takes no
`modification_parameters`: Rosemary parameterizes every atom of the protein
component, glycans and attached groups included. A config with both is
refused. For glycans on GLYCAM06j, keep the protein on ff14SB.

The audit reports glycans under Rosemary as an INFO line of
`parameters.glycan_force_field`: "glycans parameterized by Rosemary
(whole-protein mode)".

Unlike the ff14SB path, Rosemary accepts a neutral N- or C-terminus (NH2 or
COOH, with no group attached): it charges whatever chemistry it is given. A
truncated chain is still better capped with ACE or NME, so that its ends are
not charged artificially.

## Step 4: Decide on the other components

Polymers, ligands and co-solvents keep the force field they name. When they
are not on Rosemary, the audit adds an INFO line
(`parameters.rosemary_components`) that recommends Rosemary for them too.
This is not an error: Rosemary is untested in production, so you may keep
Sage for them. To follow the recommendation, name the Rosemary file on those
solutes and in `solvent.force_field`; this release applies one force field
file to every molecule but the protein.

## Step 5: Build and read the audit

```bash
polyzymd build -c config.yaml
```

The audit recomputes the NAGL charges with the declared model and checks every
term against Rosemary. Expected audit lines in `build_provenance.md`:

| Status | Check | Meaning |
|--------|-------|---------|
| WARN | `parameters.rosemary_alpha` | Rosemary is an alpha release |
| PASS | `solvent.rosemary_water` | the water is OPC3 |
| INFO or PASS | `parameters.rosemary_components` | other components on another force field, or none |

"Methods to cite" lists OPC3 (Izadi and Onufriev, 2016) and the ion set
(Sengupta et al., 2021).

PolyzyMD decides from the force field file, not its name, which charges apply:
a file that charges with NAGL takes `{method: nagl}`, a file with protein
library charges takes `{method: library}`. Every NAGL component of a build
uses one model, so the protein's model must match that of any NAGL-charged
ligand, polymer or co-solvent.
