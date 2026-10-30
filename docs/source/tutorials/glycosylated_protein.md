# From an idea and a PDB to a glycosylated protein

```{index} !glycosylation tutorial, !glycan, !modified protein
```

In this tutorial you put a glycan on a protein, build the system, check that
the build is what you asked for, and run it. You start from an idea
("chitobiose on Asn60 of ubiquitin"), a prepared PDB file and a glycan file,
and end with a simulation whose every parameter the audit has checked.

- *You learn:* how a modification is written in the config, how to see its
  plan before building, and how the build record shows what was applied.
- *You run:* ubiquitin with chitobiose (GlcNAc-β1-4-GlcNAc) N-linked to Asn60,
  in water, for 200 ps.
- *Time:* about 15 minutes. The build takes about a minute; the run takes a
  minute on a GPU and about 30 minutes on a laptop CPU.

The files are in `examples/modified_proteins`. Copy the folder and work in the
copy:

```bash
cp -r ~/polyzymd/examples/modified_proteins ~/glycan_tutorial
cd ~/glycan_tutorial
```

## Step 1: The protein

`structures/ubiquitin.pdb` is ubiquitin (PDB 1UBQ) with hydrogens added at
pH 7. PolyzyMD modifies a structure as given, so start from a protein PDB that
OpenFF can read, with every hydrogen. {doc}`prepare_pdb_for_openff` shows how
to make one from a crystal structure.

## Step 2: The glycan

`structures/glycans/glycam_G42666HT.pdb` is chitobiose, built by GLYCAM-Web
and downloaded from [GlyGen](https://www.glygen.org/) by its GlyTouCan
accession, G42666HT. Three things about this file matter:

- its residues have GLYCAM06 names (`4YB`, `0YB`), which GLYCAM06j uses to
  find their parameters;
- its first residue is `ROH`, the hydroxyl on the anomeric carbon, which
  leaves when the glycan bonds to the protein;
- it lists every bond in `CONECT` records. Files downloaded from GlyGen have
  none; {doc}`../how_to/glycosylate` shows how to add them.

## Step 3: Write the modification

Open `config_glycan.yaml`. The protein entry carries the change:

```yaml
    modifications:
      - attach:
          type: glycan
          reaction: N-glycosylation      # bonds C1 to ND2; HD22 and the ROH group leave
          chains: [A]
          glycans: {G42666HT: [N60]}     # read from structures/glycans/glycam_G42666HT.pdb
    modification_parameters:
      glycan: {force_field: [GLYCAM_06j-1.xml], charges: {method: library}}
```

`glycans:` maps an accession to the residues that get it; `N60` with
`chains: [A]` is the site `A:N60`. `modification_parameters` says which force
field each kind of change uses. Glycans use GLYCAM06j, which keeps the
ff14SB backbone of the linked Asn and gives the side chain and the sugars
their GLYCAM parameters.

## Step 4: See the plan

```bash
polyzymd validate -c config_glycan.yaml --plans
```

PolyzyMD applies the change to the PDB and prints one row per site before
anything is built:

```text
Modifications of ubq_glycan (1, checked against its PDB):
  G001  attach      A:N60      ASN->ASN  glycam_G42666HT.pdb via N-glycosylation; site G001: GLYCAM_06j-1.xml (library); +55 atoms, -1 (HD22), formal charge +0
```

`G001` is the glycan's segment. If residue 60 were not an Asn, or the
reaction did not fit, `validate` would stop here and say why.

## Step 5: Build

```bash
polyzymd build -c config_glycan.yaml
```

The build applies the change, relaxes the glycan, gives every atom its
parameters, solvates and writes `ubq_glycan-mod<hash>_300K_run1/`. The `-mod`
token in the folder name keeps a modified protein's runs apart from an
unmodified protein of the same name.

The build ends with the physics audit. Open `build_provenance.md` in the run
folder. Among its checks you will find:

```text
| PASS | `modifications.declared` | ubq_glycan | 1 changes and their parameter sets are the declared ones |
| PASS | `modifications.sites` | ubq_glycan | 62 site atoms in 1 sites, as the changes define them |
| PASS | `parameters.glycan_force_field` | ubq_glycan | 1 glycans on GLYCAM06j |
| PASS | `parameters.modified_backbone` | ubq_glycan | 6 backbone atoms of modified residues carry ff14SB charges and terms |
| PASS | `parameters.one_four` | ubq_glycan | 7073 exceptions terms match the declared sources |
```

The audit rebuilt the protein's parameters on its own and compared them with
the built system term by term. A build that fails any check stops.
{doc}`../explanation/modified_residues` explains what each of these checks
guards.

## Step 6: Run

The config runs 100 ps of NVT and 100 ps of NPT on an NVIDIA GPU
(`openmm: platform: CUDA`). On a computer without one, set `platform: "CPU"`.

```bash
polyzymd run -c config_glycan.yaml
```

The run starts from the audited build and checks that the system it runs is
the one the audit passed (`production_0/production_0_system_check.json`,
`"matches_build": true`).

## What you did

You wrote a modification as one entry of the config, saw its plan, built it,
and read the audit that proves the build carries what you declared. The same
pattern scales: the BG505 SOSIP trimer with its 63 glycans is three `attach`
entries, one per kind of chain, each with a `glycans:` map. See
{doc}`../how_to/glycosylate`.

Next:

- {doc}`../how_to/glycosylate` for more glycans, O-glycans and site tables;
- {doc}`../how_to/attach_groups` for dyes, labels and click chemistry;
- {doc}`../how_to/mutate` for mutations;
- {doc}`../how_to/analyze_modified_proteins` to select the protein and its
  attachments in analyses.
