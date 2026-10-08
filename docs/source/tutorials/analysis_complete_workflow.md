# Compare a polymer condition with water

In this tutorial you add a second condition to the quickstart study: Trp-cage
with two short chains of sulfobetaine methacrylate (SBMA). You run three
replicates of it, and then compare it with the `Water` condition in four
analyses with one command.

You learn these steps:

1. Copy a condition with `polyzymd study add-condition --from`.
2. Add polymer chains to its config.
3. Run the replicates of the new condition.
4. Run every analysis of the project with `polyzymd analyze --project`.
5. Read the comparisons with the control, and their warnings.

## Before you start

Do {doc}`first_analysis` first. This tutorial continues in its project
folder, `~/pz_quickstart`, where the `Water` condition has three replicates.

:::{admonition} Environment Setup
:class: tip

Run every command of this tutorial in the `build` environment. It holds
PACKMOL, the polymer builder and the analysis tools. From the repository
root, activate it once:

```bash
pixi shell -e build
```
:::

## Step 1: Copy the condition

```bash
cd ~/pz_quickstart
polyzymd study add-condition SBMA --from Water --study trpcage
```

The output is:

```
condition SBMA: /home/me/pz_quickstart/trpcage/conditions/sbma/config.yaml, listed in study.yaml
warning: the runs go into runs/trpcage/sbma unless you set scratch_directory in config.yaml; trajectories can use a lot of disk space, so on a cluster set scratch_directory to scratch storage
next: commit, and run polyzymd study check
```

The copy is in `trpcage/conditions/sbma/`. `Water` is the first condition of
the study, so it stays the control. Every other condition is compared with
it.

## Step 2: Add the polymer

Open `trpcage/conditions/sbma/config.yaml`. Change the name and the
description:

```yaml
name: trpcage_sbma
description: Trp-cage with two SBMA 3-mers in water with NaCl
```

Add this section at the end of the file:

```yaml
polymers:
  enabled: true
  generation_mode: dynamic
  type_prefix: SBMA
  reactions:
    initiation: default
    polymerization: default
    termination: default
  monomers:
    - label: A
      probability: 1.0
      name: SBMA
      residue_name: SBM
      smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
  length: 3
  count: 2
  packing:
    padding: 0.5
```

The section asks for two chains of three SBMA monomers. `dynamic` builds each
chain from the SMILES of the monomer, with the ATRP reaction templates that
ship with PolyzyMD. `residue_name: SBM` names the monomers in the topology.
`packing.padding: 0.5` places the chains within 0.5 nm of the protein
surface, so they touch the protein within the few picoseconds of this
tutorial. A real study uses the default, 2.0 nm, and a long simulation. For
every key, see {doc}`../how_to/dynamic_polymers`.

Check the config:

```bash
polyzymd validate -c trpcage/conditions/sbma/config.yaml
```

The summary now lists the polymer:

```
Summary:
  Name: trpcage_sbma
  Engine: openmm
  Enzyme: trpcage
  Substrate: None (apo simulation)
  Polymers: SBMA
    Count: 2
    Length: 3
    Monomer A: 100.0%
  Co-solvents: none
  Temperature: 300.0 K
  Pressure: 1.0 atm
```

The build writes the polymer fragments into `.polymer_cache/` in the current
folder. The project `.gitignore` keeps it out of git. Commit the new condition:

```bash
git add -A
git commit -m "Add the SBMA condition"
```

## Step 3: Run three replicates

```bash
polyzymd run -c trpcage/conditions/sbma/config.yaml -r 1-3
```

For each replicate, the build draws the chains, places them around the
protein with PACKMOL, and then adds water and ions. The first build also
makes the polymer fragments. The three runs take about seven minutes. The
last line is:

```
All 3 replicate(s) completed successfully.
```

Check that the study finds the runs of both conditions:

```bash
polyzymd study check trpcage
```

The condition lines are:

```
control Water: replicates [1, 2, 3] under /home/me/pz_quickstart/trpcage/conditions/water/../../../runs/trpcage/water (from config)
condition SBMA: replicates [1, 2, 3] under /home/me/pz_quickstart/trpcage/conditions/sbma/../../../runs/trpcage/sbma (from config)
```

## Step 4: Analyze the project

Open `project.yaml` and list two more analyses: the polymer contacts and the
protein-polymer hydrogen bonds.

```yaml
analyses:
  rg: {}
  rmsf: {}
  contacts: {}
  hydrogen_bonds: {}
```

Commit, and run every analysis of the project:

```bash
git add -A
git commit -m "Add the contacts and hydrogen_bonds analyses"
polyzymd analyze --project .
```

The command runs each analysis on every replicate of both conditions. It
takes less than a minute. The output is:

```
log: /home/me/pz_quickstart/logs/polyzymd-analyze-20261006-210749-pid22351.log
== study trpcage
== rg
# polyzymd analyze rg  metric mean_rg  unit A  eq 0ns  conditions 2  replicates 3,3  protocol rg/2
Water  n 3  mean 7.315  sem 0.03064  ci95 7.183 to 7.447  values 7.254, 7.346, 7.346  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
SBMA  n 3  mean 7.381  sem 0.07171  ci95 7.073 to 7.69  values 7.269, 7.515, 7.359  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
Water vs SBMA  delta +0.06592  ci95 -0.1982 to 0.33  p 0.466  p_adj 0.466  test welch_t  correction BH  family 1  d 0.6903  not_significant
warning: condition Water, condition SBMA: replicates 1, 2, 3 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected
verdict: no significant difference in mean_rg between Water and SBMA (delta +0.06592 A, 95% CI -0.1982 to 0.33, p_adj 0.466, p 0.466, n 3 vs 3)

== rmsf
# polyzymd analyze rmsf  metric core_rmsf  unit A  run core_rmsf  eq 0ns  conditions 2  replicates 3,3  protocol rmsf/2
Water  n 3  mean 0.3463  sem 0.01661  ci95 0.2748 to 0.4178  values 0.3669, 0.3585, 0.3134
SBMA  n 3  mean 0.3687  sem 0.03553  ci95 0.2159 to 0.5216  values 0.3071, 0.4302, 0.369
Water vs SBMA  delta +0.02248  ci95 -0.1066 to 0.1515  p 0.6089  p_adj 0.6089  test welch_t  correction BH  family 1  d 0.4679  not_significant
verdict: no significant difference in core_rmsf between Water and SBMA (delta +0.02248 A, 95% CI -0.1066 to 0.1515, p_adj 0.6089, p 0.6089, n 3 vs 3)

== contacts
# polyzymd analyze contacts  metric coverage  unit none  run coverage  eq 0ns  conditions 2  replicates 3,3  protocol contacts/2
Water  n 3  mean 0  sem 0  ci95 na  values 0, 0, 0
SBMA  n 3  mean 0.1  sem 0.02887  ci95 -0.02421 to 0.2242  values 0.1, 0.15, 0.05
Water vs SBMA  delta +0.1  ci95 -0.02421 to 0.2242  p na  p_adj na  test welch_t  correction BH  d 2.828  not_testable
warning: condition Water has the same coverage in every replicate, so its interval is not estimable
warning: the 95 percent interval of condition SBMA extends past the bounds 0 to 1 of coverage, where a t interval is not reliable
warning: contacts: polymer_selection 'chainid C' matched no atoms in Water replicate 1, 2, 3, so contact there is 0 (none of those atoms to touch) and no comparison with it is tested. Check the selection if that condition has them.
verdict: not testable: coverage for Water vs SBMA, as Water has no partner: polymer_selection 'chainid C' matched no atoms in replicate 1, 2, 3 (n 3 vs 3)

== hydrogen_bonds
# polyzymd analyze hydrogen_bonds  metric protein_polymer_mean_hbonds  unit none  run protein_polymer_mean_hbonds  eq 0ns  conditions 2  replicates 3,3  protocol hydrogen_bonds/2
Water  n 3  mean 0  sem 0  ci95 na  values 0, 0, 0
SBMA  n 3  mean 1.083  sem 0.5069  ci95 -1.098 to 3.264  values 2, 0.25, 1
Water vs SBMA  delta +1.083  ci95 -1.098 to 3.264  p na  p_adj na  test welch_t  correction BH  d 1.745  not_testable
warning: condition Water has the same protein_polymer_mean_hbonds in every replicate, so its interval is not estimable
warning: the 95 percent interval of condition SBMA extends past the bounds 0 to inf of protein_polymer_mean_hbonds, where a t interval is not reliable
warning: hydrogen_bonds: second group 'chainid C' matched no atoms in Water replicate 1, 2, 3, so the hydrogen-bond count there is 0 (none of those atoms to touch) and no comparison with it is tested. Check the selection if that condition has them.
verdict: not testable: protein_polymer_mean_hbonds for Water vs SBMA, as Water has no partner: second group 'chainid C' matched no atoms in replicate 1, 2, 3 (n 3 vs 3)
```

Your values differ, because each run adds up the forces in a different
order. The outputs on these pages come from separate runs, so `Water`
replicate 1 here does not match the quickstart output exactly.

These runs are three replicates of four frames. A difference must be large to
be significant with so few samples, and a result that is not significant does
not show that the conditions are the same. The numbers show the steps, not a
result about SBMA.

## Step 5: Read the comparisons

Each analysis prints one line per condition, then one comparison line,
`Water vs SBMA`. The comparison gives the difference of the means (`delta`)
with its 95 % interval, the p value of Welch's t test, the p value after the
{term}`Benjamini-Hochberg` correction (`p_adj`) and the effect size `d`. The
last word says whether the difference is significant.

- **rg and rmsf.** The radius of gyration and the fluctuation of Trp-cage
  show no significant difference here.
- **contacts.** `coverage` is the fraction of protein residues that the
  polymer touches on at least one frame. In this run, the chains in `SBMA`
  touch a few percent of the residues in each replicate. `Water` has no
  polymer, so its coverage is 0.
- **hydrogen_bonds.** In this run, the chains in `SBMA` form about one
  hydrogen bond with the protein per frame. In a run this short, one
  replicate can show 0.

Read every `warning:` line. Each says what limits the result:

- `matched no atoms in Water` says that the control has no chain C, so its
  value is 0 by definition. That 0 is not a measurement, so the comparison
  with `Water` is `not testable`. The warning asks you to check that this is
  expected. Here it is.
- `the same ... in every replicate` says that `Water` has no variance. The
  test then has little power, and the verdict says so.
- `extends past the bounds` says that a t interval does not suit the values,
  such as a fraction close to 0.
- `fewer than 20 effective samples` comes from the four frames of each
  replicate. A real study has many more.

## Step 6: Find the results

Each analysis has its folder in `trpcage/results/`:

```text
trpcage/results/
├── rg/
├── rmsf/
├── contacts/
│   ├── report.json
│   ├── figures/contacts/
│   │   ├── contacts_coverage_comparison.png
│   │   ├── contacts_contact_fraction_profile.png
│   │   ├── contacts_contact_fraction_difference.png
│   │   └── contacts_class_bars.png
│   └── polyzymd_results/residue_occlusion/
└── hydrogen_bonds/
    ├── report.json
    ├── figures/hydrogen_bonds/
    │   └── hbonds_protein_polymer_mean_hbonds_comparison.png
    └── polyzymd_results/hydrogen_bonds_protein_polymer/
```

Each `polyzymd_results/` folder holds one folder per condition, with one
folder per replicate. Each replicate folder holds the values and the
`record.json` that says how they were made.

To compare the residues one by one, report another result of contacts:

```bash
polyzymd analyze contacts --study trpcage --run contact_fraction_residues
```

It reuses the stored values. Its comparison line is:

```
Water vs SBMA  labels 20  tested 5  family 5  test welch_t  correction BH  lower 0  higher 0
```

Of the 20 residues, 5 have values that vary, so 5 are tested. None differs
after the correction. A `warning:` line follows for each residue whose test
or interval is limited. The JSON report (`--format json`) holds every row.

## What you did

You added a polymer condition to a study, ran it, and compared it with the
control in four analyses with one command. Next, measure how much of the
protein surface the polymer covers in {doc}`sasa_analysis`. For the settings
of each analysis, see {doc}`../how_to/analysis_contacts_quickstart` and
{doc}`../how_to/hydrogen_bonds`. For every line of the report, see
{doc}`../how_to/analysis_agent_protocol`.
