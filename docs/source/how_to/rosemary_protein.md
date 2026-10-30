# Put a protein on OpenFF Rosemary alpha

This guide sets up a protein on OpenFF Rosemary alpha in OPC3 water. Rosemary
is OpenFF's protein force field in development; the alpha is not for
production. The fields are listed in
[OpenFF Rosemary alpha](protein-rosemary-alpha) of the configuration
reference.

:::{note}
Run `polyzymd` from the `build` pixi environment (`pixi run -e build ...`).
:::

## 1. Name Rosemary on the protein

```yaml
solutes:
  - type: protein
    pdb_path: "structures/enzyme.pdb"
    force_field: ["openff_no_water-3.0.0-alpha2b.offxml"]
    charges: {method: nagl, model: openff-gnn-am1bcc-1.0.0.pt}
```

Rosemary parameterizes the whole protein component, mutations and
modifications included, and NAGL charges the whole molecule.

## 2. Use OPC3 water

OpenFF pairs Rosemary with OPC3 water. Set the water model:

```yaml
solvent:
  force_field: ["openff-2.3.0.offxml"]
  water:
    model: "opc3"
```

PolyzyMD loads the force fields in this order: the protein's (Rosemary), then
`solvent.force_field`, then `opc3.offxml`. The water and the Na+ and Cl- ions
take their parameters from `opc3.offxml` (OPC3 water and the Li/Merz 12-6 ions
fitted for it); see [Water Models](water-models).

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

## 3. Decide on the other components

Polymers, ligands and co-solvents keep the force field they name. When they
are not on Rosemary, the audit adds an INFO line
(`parameters.rosemary_components`) that recommends Rosemary for them too.
This is not an error: Rosemary is untested in production, so you may keep
Sage for them. To follow the recommendation, name the Rosemary file on those
solutes and in `solvent.force_field`; this release applies one force field
file to every molecule but the protein.

Glycans attached to the protein are part of the protein component, so
Rosemary parameterizes them too. The audit reports this as an INFO line of
`parameters.glycan_force_field`: "glycans parameterized by Rosemary
(whole-protein mode)". A glycosylated protein on Rosemary takes no
`modification_parameters`; for glycans on GLYCAM06j, keep the protein on
ff14SB.

## 4. Build and read the audit

```bash
polyzymd build -c config.yaml
```

Expected audit lines in `build_provenance.md`:

| Status | Check | Meaning |
|--------|-------|---------|
| WARN | `parameters.rosemary_alpha` | Rosemary is an alpha release |
| PASS | `solvent.rosemary_water` | the water is OPC3 |
| INFO or PASS | `parameters.rosemary_components` | other components on another force field, or none |

"Methods to cite" lists OPC3 (Izadi and Onufriev, 2016) and the ion set
(Sengupta et al., 2021).
