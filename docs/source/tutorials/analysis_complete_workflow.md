# Compare a polymer condition with water

:::{admonition} At a glance
:class: note

- **You will learn:** how a study compares each condition with its control,
  and how to read the comparison and the warnings that limit it.
- **You need:** the project of {doc}`first_analysis`, with three replicates
  of `Water`.
- **Time:** about 20 minutes on a laptop CPU, 5 to 7 of them building and
  simulating.
:::

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

:::{admonition} Environment setup
:class: tip

Run every command of this tutorial in the `build` environment. It holds
PACKMOL, the polymer builder and the analysis tools. From the repository
root, activate it once:

```bash
pixi shell -e build
```

**On the NVIDIA GPU of this machine:** in step 3, run `polyzymd build` with
the same `-c` and `-r` in `build`, then `polyzymd run` in your `sim-cuda-*`
shell. The copy keeps the `openmm:` settings of `Water`.

**With GROMACS:** if `Water` uses `config_gromacs.yaml`, the copy keeps
`engine: gromacs`, and every step works the same. `validate` then prints
`Engine: gromacs`.
:::

## Step 1: Copy the condition

```bash
cd ~/pz_quickstart
polyzymd study add-condition SBMA --from Water --study trpcage
```

The output is:

```
condition SBMA: /home/me/pz_quickstart/trpcage/conditions/sbma/config.yaml, listed in study.yaml
name: set to 'sbma' in the copy; change description: and the settings in which SBMA differs from Water
warning: the runs go into runs/trpcage/sbma unless you set scratch_directory in config.yaml; trajectories can use a lot of disk space, so on a cluster set scratch_directory to scratch storage
next: commit, and run polyzymd study check
```

The copy is in `trpcage/conditions/sbma/`. `Water` is the first condition of
the study, so it stays the control. Every other condition is compared with
it.

## Step 2: Add the polymer

Open `trpcage/conditions/sbma/config.yaml`. Its `name:` is `sbma`. Change
it to `trpcage_sbma`, and change the description:

```yaml
name: trpcage_sbma
description: Trp-cage with two SBMA 5-mers in water with NaCl
```

Add the polymer entry below to the `solutes:` list, after the protein entry.
The protein entry is shown for context; do not copy it a second time.

:::{tip}
Indent `- type: polymer` exactly like `- type: protein` above it, and the
lines below it by the same amount. Otherwise `validate` stops with
`not valid YAML: expected <block end>, but found '-'`.
:::

```yaml
solutes:
  - type: protein
    name: "trpcage"
    # ... the rest of the protein entry, unchanged
  - type: polymer
    generation_mode: dynamic
    type_prefix: SBMA
    monomers:
      - label: A
        probability: 1.0
        name: SBMA
        residue_name: SBM
        smiles: "[*]CC([*])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"
    end_groups:
      head: [{label: A, name: SBMA, residue_name: SBM, smiles: "[*]C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]", probability: 1.0}]
      tail: [{label: A, name: SBMA, residue_name: SBM, smiles: "[*]C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]", probability: 1.0}]
    length: 5
    count: 2
    force_field: ["openff-2.3.0.offxml"]
    charges: {method: nagl, model: "openff-gnn-am1bcc-1.0.0.pt"}
    packing:
      padding: 0.5
```

The entry asks for two chains of five SBMA monomers. `dynamic` builds each
chain with mBuild from the repeat unit of the monomer: its SMILES has `*` on
the two carbons that bond to the neighbouring units. `end_groups` puts an
SBMA unit that kept its C=C double bond at each end, as in the chains PolyzyMD built before 1.3.0,
so a 5-mer is an end, three repeat units and an end. `residue_name: SBM` names the
monomers in the topology.
`force_field` is the same file as `solvent.force_field`. `charges` charges
each new chain with {term}`NAGL` and the model you name.
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
    Length: 5
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

```{figure} ../_static/renders/polymer_system.webp
:figclass: pz-render
:class: dark-light
:alt: Two coiled polymer chains drawn as thick blue sticks with red, blue and yellow sulfonate and ammonium groups, beside a small protein drawn as a gray cartoon, in water shown as small light spheres with a few purple sodium and yellow chloride ions.

What you just built: the two SBMA 5-mers (carbons blue; the sulfonate S
yellow and O red, the ammonium N blue) beside Trp-cage (gray), in water with
Na⁺ (purple) and Cl⁻ (yellow). The coordinates are those of the build, before
the simulation; water and ions in front of the chains and the protein are cut
away. The shipped example `examples/polymer` builds the same system.
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
protein-polymer hydrogen bonds. Both find the polymer with the selection
`chainid C`, because the build puts the protein in chain A, a small molecule
in chain B, every polymer chain in chain C and the solvent in chains D and
later.

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
log: /home/me/pz_quickstart/logs/polyzymd-analyze-20261026-113501-pid274155.log
== study trpcage
== rg
note: replaced the stored results/rg/report.json, which reported mean_rg
# polyzymd analyze rg  metric mean_rg  unit A  eq 0ns  conditions 2  replicates 3,3  protocol rg/2
Water  n 3  mean 7.294  sem 0.0631  ci95 7.023 to 7.566  values 7.391, 7.316, 7.176  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.001 ns
SBMA  n 3  mean 7.277  sem 0.04271  ci95 7.093 to 7.461  values 7.243, 7.362, 7.227  replicates 1, 2, 3  g 1, 1, 1  n_eff 4, 4, 4  eq_detected 0.002 ns
Water vs SBMA  delta -0.01697  ci95 -0.2405 to 0.2066  p 0.8361  p_adj 0.8361  test welch_t  correction BH  family 1  d -0.1818  not_significant
warning: condition Water, condition SBMA: replicates 1, 2, 3 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected
verdict: no significant difference in mean_rg between Water and SBMA (delta -0.01697 A, 95% CI -0.2405 to 0.2066, p_adj 0.8361, p 0.8361, n 3 vs 3)

== rmsf
note: replaced the stored results/rmsf/report.json, which reported core_rmsf
# polyzymd analyze rmsf  metric core_rmsf  unit A  run core_rmsf  eq 0ns  conditions 2  replicates 3,3  protocol rmsf/2
Water  n 3  mean 0.3375  sem 0.02121  ci95 0.2462 to 0.4287  values 0.3571, 0.3602, 0.2951
SBMA  n 3  mean 0.3158  sem 0.01206  ci95 0.2639 to 0.3677  values 0.3392, 0.2991, 0.3091
Water vs SBMA  delta -0.02165  ci95 -0.097 to 0.0537  p 0.4371  p_adj 0.4371  test welch_t  correction BH  family 1  d -0.7243  not_significant
verdict: no significant difference in core_rmsf between Water and SBMA (delta -0.02165 A, 95% CI -0.097 to 0.0537, p_adj 0.4371, p 0.4371, n 3 vs 3)

== contacts
# polyzymd analyze contacts  metric coverage  unit none  run coverage  eq 0ns  conditions 2  replicates 3,3  protocol contacts/2
Water  n 3  mean 0  sem 0  ci95 na  values 0, 0, 0
SBMA  n 3  mean 0.1333  sem 0.0441  ci95 -0.0564 to 0.3231  values 0.05, 0.2, 0.15
Water vs SBMA  delta +0.1333  ci95 -0.0564 to 0.3231  p na  p_adj na  test welch_t  correction BH  d 2.469  not_testable
warning: condition Water has the same coverage in every replicate, so its interval is not estimable
warning: the 95 percent interval of condition SBMA extends past the bounds 0 to 1 of coverage, where a t interval is not reliable
warning: contacts: polymer_selection 'chainid C' matched no atoms in Water replicate 1, 2, 3, so contact there is 0 (none of those atoms to touch) and no comparison with it is tested. Check the selection if that condition has them.
verdict: not testable: coverage for Water vs SBMA, as Water has no partner: polymer_selection 'chainid C' matched no atoms in replicate 1, 2, 3 (n 3 vs 3)

== hydrogen_bonds
# polyzymd analyze hydrogen_bonds  metric protein_polymer_mean_hbonds  unit none  run protein_polymer_mean_hbonds  eq 0ns  conditions 2  replicates 3,3  protocol hydrogen_bonds/2
Water  n 3  mean 0  sem 0  ci95 na  values 0, 0, 0
SBMA  n 3  mean 0.4167  sem 0.2205  ci95 -0.532 to 1.365  values 0.5, 0.75, 0
Water vs SBMA  delta +0.4167  ci95 -0.532 to 1.365  p na  p_adj na  test welch_t  correction BH  d 1.543  not_testable
warning: condition Water has the same protein_polymer_mean_hbonds in every replicate, so its interval is not estimable
warning: the 95 percent interval of condition SBMA extends past the bounds 0 to inf of protein_polymer_mean_hbonds, where a t interval is not reliable
warning: hydrogen_bonds: second group 'chainid C' matched no atoms in Water replicate 1, 2, 3, so the hydrogen-bond count there is 0 (none of those atoms to touch) and no comparison with it is tested. Check the selection if that condition has them.
verdict: not testable: protein_polymer_mean_hbonds for Water vs SBMA, as Water has no partner: second group 'chainid C' matched no atoms in replicate 1, 2, 3 (n 3 vs 3)
```

Your values differ, because each run adds up the forces in a different
order. The outputs on these pages come from separate runs, so the `Water`
values here do not match the earlier tutorials exactly.

The two `note:` lines say that the stored reports of `rg` and `rmsf`, which
had only `Water`, are replaced by reports with both conditions. The values of
`Water` are not measured again: the analyses reuse the values stored by the
earlier tutorials.

These runs are three replicates of four frames. A difference must be large to
be significant with so few samples, and a result that is not significant does
not show that the conditions are the same. The numbers show the steps, not a
result about SBMA.

## Step 5: Read the comparisons

Each analysis prints one line per condition, then one comparison line,
`Water vs SBMA`. The comparison gives the difference of the means (`delta`)
with its 95 % interval, the p value of Welch's t test, the p value after the
{term}`Benjamini-Hochberg` correction (`p_adj`) and the effect size `d`. The
last word says whether the difference is significant, or `not_testable` when
no test can be made.

- **rg and rmsf.** The radius of gyration and the fluctuation of Trp-cage
  show no significant difference here.
- **contacts.** `coverage` is the fraction of protein residues that the
  polymer touches on at least one frame. In this run, the chains in `SBMA`
  touch 1 to 4 of the 20 residues (5 to 20 %) in each replicate. `Water` has
  no polymer, so its coverage is 0 and the comparison is `not_testable`.
- **hydrogen_bonds.** In this run, the chains in `SBMA` form fewer than one
  hydrogen bond with the protein per frame, and replicate 3 forms none. In a
  run this short, one replicate can show 0. As for contacts, `Water` has no
  polymer, so the comparison is `not_testable`.

Read every `warning:` line. Each says what limits the result:

- `matched no atoms in Water` says that the control has no chain C, so its
  value is 0 by definition. That 0 is not a measurement, so the comparison
  with `Water` is `not testable`. The warning asks you to check that this is
  expected. Here it is.
- `the same ... in every replicate` says that `Water` has no variance, so it
  has no interval.
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

It reuses the stored values. Its comparison line and its verdict are:

```
Water vs SBMA  labels 20  tested 0  family na  test welch_t  correction BH  lower 0  higher 0
...
verdict: not testable: contact_fraction for Water vs SBMA, as Water has no partner: polymer_selection 'chainid C' matched no atoms in replicate 1, 2, 3 (n 3 vs 3)
```

None of the 20 residues is tested: as for `coverage`, `Water` has no polymer,
so its contact fraction is 0 by definition. This command also writes the per-residue figures to
`figures/contacts/`: `contacts_contact_fraction_profile.png` and, because the
study has two conditions, `contacts_contact_fraction_difference.png`. The
default result, `coverage`, writes only the two figures above. A `warning:` line follows for each residue whose test
or interval is limited. The JSON report (`--format json`) holds every row.

## What you did

You added a polymer condition to a study, ran it, and compared it with the
control in four analyses with one command. Next, measure how much of the
protein surface the polymer covers in {doc}`sasa_analysis`. For the settings
of each analysis, see {doc}`../how_to/analysis_contacts_quickstart` and
{doc}`../how_to/hydrogen_bonds`. For every line of the report, see
{doc}`../how_to/analysis_agent_protocol`.
