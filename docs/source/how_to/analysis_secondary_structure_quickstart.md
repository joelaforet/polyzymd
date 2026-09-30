# Secondary structure analysis: quick start

Assign each residue a DSSP secondary-structure class on every production frame
of every replicate, and compare how much of the protein is helix, strand or
coil, overall or residue by residue, with the replicate as the sampling unit.

```{versionadded} 1.3.0
Secondary structure analysis was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the measurement?** For what each shipped function
measures, see {doc}`../reference/analysis_functions`; for the statistics, see
{doc}`../explanation/analysis_statistics_best_practices`.
```

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

## What is measured

On every production frame, `mdtraj.compute_dssp(simplified=False)` assigns
each residue of the selection one of the eight DSSP classes of Kabsch and
Sander (1983) from its backbone hydrogen bonds and geometry:

| Class | DSSP code | Group |
|---|---|---|
| `alpha_helix` | H | `helix` |
| `3_10_helix` | G | `helix` |
| `pi_helix` | I | `helix` |
| `extended_strand` | E | `strand` |
| `isolated_bridge` | B | `strand` |
| `turn` | T | `coil` |
| `bend` | S | `coil` |
| `loop` | blank | `coil` |
| `unassigned` | NA | none |

The groups are the three classes of MDTraj's simplified DSSP. MDTraj gives
`NA` to a residue it cannot assign, for example one missing a backbone atom;
PolyzyMD counts those as `unassigned`, in no group, and warns. For each
replicate, each residue's value of a class or group is the fraction of
production frames it spends in it.

The selection must hold whole residues, because DSSP reads every backbone atom.
Each chain of the topology stays its own chain, so DSSP never pairs residues of
different chains. The coordinates are used as loaded; make a protein split
across a periodic boundary whole before relying on the assignment.

## From the command line

```bash
polyzymd analyze secondary_structure -c noPoly/config.yaml -c SBMA50/config.yaml \
  --label "No polymer" --label "SBMA 50%" --eq 200ns
```

The first `-c` is the control. One pass over each replicate gives every class
and group. By default the report shows `helix`: for each replicate, the
fraction of residue-frames in the helix group, which is the mean over residues
of each residue's helix fraction. The replicate values are summarised per
condition, and every other condition is compared with the control by Welch's t
test with the Benjamini-Hochberg correction. Pick another result with `--run`:

| `--run` | One value per replicate |
|---|---|
| `helix`, `strand`, `coil` (default `helix`) | Fraction of residue-frames in the group |
| any class, such as `alpha_helix` or `unassigned` | Fraction of residue-frames in the class |
| `<name>_residues` | Each residue's fraction of frames in the class or group, compared residue by residue |

A per-residue comparison is corrected over every residue of every compared
condition, and the text report gives, for each condition, how many residues
are significantly lower and higher than in the control and lists them. Every
per-residue row is kept in the JSON report.

The one setting, passed with `--set`:

| Setting | Default | Meaning |
|---|---|---|
| `selection` | `protein` | Whole residues to assign |

For example, to compare the helix content of one domain, residue by residue:

```bash
polyzymd analyze secondary_structure -c noPoly/config.yaml -c SBMA50/config.yaml \
  --eq 200ns --set "selection=protein and resid 1:120" --run helix_residues
```

Add `--stride 5` to assign every fifth production frame, `--format json` for
the full report, `--replicates 1-3` to use only some replicates, and
`--recompute` to ignore stored results. The selection is recorded under
`provenance.settings` in the JSON report.

```{note}
The plugin used before this version reported only the three groups, counted
unassigned residues as coil without a warning, merged every chain into one,
and defaulted to `protein and chainid A`. Its helix, strand and coil fractions
equal this version's for a single-chain protein with no unassigned residues.
```

## Figures

`polyzymd analyze secondary_structure` writes these figures to
`<output-dir>/figures/secondary_structure/`; `--no-plots` skips them.

| Figure | What it shows |
|---|---|
| `ss_content_bars` | The helix, strand and coil fractions of every condition, with every replicate value |
| `ss_<name>_comparison` | For a total: each condition's mean with its interval and every replicate value |
| `ss_<name>_profile` | For `<name>_residues`: each residue's fraction per replicate and each condition's mean with its interval |
| `ss_groups_<name>` | For `<name>_residues`: one panel per condition with each residue's helix, strand and coil fractions |
| `ss_<name>_difference` | For `<name>_residues` with several conditions: each condition minus the control at every residue, with the interval of the difference and the significant residues marked |

## From Python

```python
import polyzymd as pz
from polyzymd.analyses.functions import DSSP_PARTS, dssp_occupancy

study = pz.Study.from_configs(
    {"No polymer": "noPoly/config.yaml", "SBMA 50%": "SBMA50/config.yaml"},
    equilibration="200ns",
)
rows = study.per_replicate(
    dssp_occupancy,
    pz.select("protein"),
    unit=None,
    labels=lambda u: u.select_atoms("protein").residues.resids,
    parts=DSSP_PARTS,
    bounds=(0.0, 1.0),
)
helix = rows["helix"].over_labels("mean", "helix_fraction")
print(helix.compare(control="No polymer").to_agent_text())
print(rows["helix"].compare(control="No polymer").to_agent_text())  # residue by residue
```

`dssp_occupancy(atoms, frames)` returns one row per name in `DSSP_PARTS`, the
eight classes and `unassigned` from `DSSP_CLASSES`, then the groups of
`DSSP_GROUPS`, with one column per residue. `parts=` turns each row into its own
result, and `over_labels("mean")` turns a residue profile into one value per
replicate.

## References

**Kabsch W, Sander C.** (1983) "Dictionary of protein secondary structure:
pattern recognition of hydrogen-bonded and geometrical features."
*Biopolymers* 22:2577-2637. https://doi.org/10.1002/bip.360221211

## Next steps

- **What each function measures**: {doc}`../reference/analysis_functions`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **RMSF analysis**: {doc}`analysis_rmsf_quickstart`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
