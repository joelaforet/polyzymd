# Build a Protein Conjugate

This tutorial walks through the current PolyzyMD conjugation workflow from the
user's point of view: start with a cleaned, unmodified protein PDB, describe the
covalent attachment you want, run the public Python API, and inspect the files
PolyzyMD writes.

By the end, you will know how to:

- describe an attachment as **site + moiety + mechanism**;
- choose between the built-in `nhs_lys` and `n_glycosylation` mechanisms;
- call `build_conjugate()`, `build_conjugate_from_config()`, or
  `ConjugationEngine` from public API code;
- read the main `ConjugationResult` artifact paths; and
- reason about custom SMARTS without assuming unsupported plugin behavior.

```{important}
Run conjugation examples inside a PolyzyMD pixi environment. The workflow uses
heavy simulation and chemistry packages such as OpenMM, OpenFF, Pablo, and
RDKit-backed tooling.
```

## Before You Start

You need:

- a cleaned **unmodified** protein PDB, for example
  `structures/protein_clean.pdb`;
- a residue number and chain ID for the attachment site;
- either a polymer recipe for an NHS-lysine polymer conjugate, or a SMILES string
  for an N-glycosylation moiety; and
- an output directory where PolyzyMD can write construction artifacts.

If you are starting from a raw crystal structure, first work through
{doc}`prepare_pdb_for_openff`. For field-level lookup while you are learning,
see {doc}`../reference/protein_modification_config`.

The important modeling assumption in this tutorial is that the protein PDB is
the **source state**. Do not pre-edit the lysine or asparagine into a product
residue yourself. PolyzyMD needs the unmodified residue so the mechanism can
find the target atom and the leaving atoms.

## Step 1: Describe One Attachment

Every conjugation request has the same conceptual shape:

```text
clean protein PDB
  + site       -> where on the protein to attach
  + moiety     -> what to attach
  + mechanism  -> how the bond is made
  -> conjugated product artifacts
```

For example, an N-glycosylation attachment to asparagine 67 can be written as a
plain Python dictionary and validated by `ConjugateBuildRequest`:

```python
attachment = {
    "name": "asn67-glcnac",
    "site": {
        "chain_id": "A",
        "residue_name": "ASN",
        "residue_number": 67,
        # atom_name is optional here: n_glycosylation knows the target is ND2.
    },
    "moiety": {
        "name": "glcnac",
        "smiles": "OC1C(O)C(O)C(CO)OC1NC(C)=O",
        "residue_name": "NAG",
    },
    "mechanism": {"name": "n_glycosylation"},
}
```

For SMILES-defined moieties, `residue_name` is required. PolyzyMD treats the
entire SMILES moiety as **one residue** in the generated product. This is useful
for a small glycan-like group or label. Polymers are different: a polymer recipe
generates a multi-residue moiety, usually with one reactive monomer selected for
the covalent bond.

## Step 2: Run the Public API for a SMILES Moiety

Create `scripts/build_asn67_glcnac.py`:

```python
from pathlib import Path

from polyzymd.builders.conjugation import ConjugateBuildRequest, build_conjugate


request = ConjugateBuildRequest(
    protein_pdb_path=Path("structures/protein_clean.pdb"),
    output_dir=Path("artifacts/asn67-glcnac"),
    free_polymer_seed=17,
    attachments=(
        {
            "name": "asn67-glcnac",
            "site": {
                "chain_id": "A",
                "residue_name": "ASN",
                "residue_number": 67,
            },
            "moiety": {
                "name": "glcnac",
                "smiles": "OC1C(O)C(O)C(CO)OC1NC(C)=O",
                "residue_name": "NAG",
            },
            "mechanism": {"name": "n_glycosylation"},
        },
    ),
)

result = build_conjugate(request)

print(result.status)
print(result.crosslinked_conjugate_pdb_path)
print(result.minimized_conjugate_pdb_path)
print(result.solvated_pdb_path)
```

Run it with pixi:

```bash
pixi run -e build python scripts/build_asn67_glcnac.py
```

You can make the engine object explicit if you prefer an object-oriented style:

```python
from polyzymd.builders.conjugation import ConjugationEngine

engine = ConjugationEngine()
result = engine.build(request)
```

Both forms go through the same public orchestration boundary.

## Step 3: Inspect the Result

`ConjugationResult` is a lightweight summary of the construction workflow. The
main fields are paths to files you can open in a viewer or pass to later build
steps:

```python
for label, path in result.artifact_paths.items():
    print(f"{label}: {path}")
```

Common paths include:

| Field | What to inspect |
|-------|-----------------|
| `crosslinked_conjugate_pdb_path` | Product-state PDB after covalent assembly. |
| `minimized_conjugate_pdb_path` | OpenMM-relaxed product from the construction smoke workflow, when written. |
| `relaxed_conjugate_pdb_path` | Relaxed conjugate path from workflows that expose a separate relaxed product. |
| `solvated_pdb_path` | Solvated system PDB after the conjugate is placed in a box. |
| `workflow_json_path` | JSON sidecar describing the workflow and assumptions. |

Some fields can be `None` depending on the mechanism and settings. A useful
first success check is that the crosslinked and minimized PDB paths exist, and
that the modified residue names and new linkage are visible in a molecular
viewer.

Final relaxation/minimization is part of PolyzyMD's OpenMM/Pablo product-state
workflow. RDKit may help create or inspect molecular inputs, but there is no
RDKit final-minimization fallback for the conjugated product.

## Step 4: Build an NHS-Lys Polymer Conjugate from Config

Use the NHS-lysine workflow when you have an NHS-activated polymer and want to
couple it to a lysine side-chain amine.

In this path, the protein still starts unmodified. The attachment tells PolyzyMD
which lysine to target and gives a polymer recipe with one NHS-containing
reactive monomer:

```yaml
conjugation:
  enabled: true
  mode: construct
  attachments:
    - name: lys23-sbma-egpma-nhs
      site:
        chain_id: A
        residue_name: LYS
        residue_number: 23
        # atom_name can be omitted for the built-in NHS-Lys defaults.
      moiety:
        name: SBMA-EGPMA-NHS
        recipe:
          name: SBMA-EGPMA-NHS
          length: 9
          seed: 7
          forced_reactive_monomer_label: C
          monomers:
            - label: A
              name: SBMA
              residue_name: SBM
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
              probability: 0.945
            - label: B
              name: EGPMA
              residue_name: EGP
              smiles: "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
              probability: 0.045
            - label: C
              name: NHS
              residue_name: NHS
              smiles: "CC(=C)C(=O)ON1C(=O)CCC1=O"
              probability: 0.010
      mechanism:
        name: nhs_lys_amide
```

This snippet is only the `conjugation:` block. In a real `config.yaml`, keep the
normal `name`, `engine`, `enzyme`, `solvent`, `thermodynamics`, and simulation
phase sections as you would for any PolyzyMD build.

Call the config-driven public API:

```python
from pathlib import Path

from polyzymd.builders.conjugation import build_conjugate_from_config


result = build_conjugate_from_config(
    Path("config.yaml"),
    output_dir=Path("artifacts/lys23-polymer"),
    free_polymer_seed=7,
)

print(result.artifact_paths)
```

Run the script with:

```bash
pixi run -e build python scripts/build_lys23_polymer.py
```

## Step 5: Choose the Built-In Mechanism

For your first conjugation runs, choose one of the built-in mechanisms instead
of writing custom chemistry.

| Use this mechanism | When your request looks like | Current public UX |
|--------------------|------------------------------|-------------------|
| `nhs_lys` or `nhs_lys_amide` | An NHS ester on a polymer reacts with lysine `NZ` to form an amide. | Config-driven polymer recipe with an NHS reactive monomer. |
| `n_glycosylation` | A reducing-end glycan-like SMILES moiety reacts with asparagine `ND2`. | Direct `ConjugateBuildRequest` from a cleaned protein plus one or more SMILES moieties. |

Use `nhs_lys` when the thing you are adding is a polymer made of multiple
monomer residues. Use `n_glycosylation` when the thing you are adding can be
represented as a single SMILES-defined moiety residue, such as a small glycan
fragment in the current workflow.

In both cases, the mechanism can fill in known target atoms. For example,
`n_glycosylation` targets asparagine `ND2`, so you can omit
`site.atom_name` in the attachment. Add `atom_name` only when you need a strict
site check or are working with an advanced explicit configuration.

## Step 6: Read SMARTS as a Map of Roles

Once the built-ins make sense, SMARTS become easier to reason about. A reaction
SMARTS describes which atoms are consumed and produced. PolyzyMD's mechanism
metadata adds user-facing roles to the atom-map numbers:

```yaml
mechanism:
  name: example_amide_metadata
  reaction_smarts: "[N:1]([H:2]).[C:3](=[O:4])[O:5]>>[N:1][C:3](=[O:4])"
  atom_roles:
    - map_number: 1
      participant: site
      role: linking
      label: protein_nucleophile
    - map_number: 2
      participant: site
      role: leaving
      label: site_hydrogen
    - map_number: 3
      participant: moiety
      role: linking
      label: moiety_acyl_carbon
    - map_number: 4
      participant: moiety
      role: geometry_anchor
      label: carbonyl_oxygen
    - map_number: 5
      participant: moiety
      role: leaving
      label: moiety_leaving_oxygen
```

Read this as:

1. atom map `1` is the protein atom that will form the new bond;
2. atom map `3` is the moiety atom that will form the other side of the bond;
3. atom maps `2` and `5` are removed during bond formation; and
4. atom map `4` is retained and can help orient or validate the local chemistry.

The built-in templates use the same idea. For NHS-Lys, the linking atoms are the
lysine nitrogen and NHS acyl carbon. For N-glycosylation, the linking atoms are
asparagine `ND2` and the glycan anomeric carbon; the mechanism also identifies
the anomeric hydroxyl atoms that leave.

## Step 7: Validate Custom SMARTS Conservatively

Custom SMARTS are currently best treated as **mechanism metadata and validation
input**, not as an arbitrary external plugin system. The executable public paths
are still centered on the built-in mechanisms above. Supplying a new
`reaction_smarts` string does not, by itself, teach PolyzyMD how to place the
moiety, remove atoms from a product PDB, parameterize a new residue, or run a
new coordinate-surgery backend.

When you design a future custom mechanism, use this checklist:

- Start with the smallest realistic reactants, not the full protein.
- Atom-map every atom whose bonding changes.
- Mark exactly two `linking` atoms: one on `site`, one on `moiety`.
- Mark atoms that disappear as `leaving` and keep them participant-scoped.
- Include at least one retained `geometry_anchor` near the reactive center when
  orientation matters.
- Confirm that the target protein residue has the expected atom and removable
  hydrogens in the cleaned PDB.
- Decide the product residue names and force-field/template strategy before
  running expensive builds.

If your chemistry is not NHS-Lys or the current N-glycosylation path, expect to
write or request a supported reaction template rather than relying on SMARTS
alone.

## What You Learned

You have followed the current conjugation workflow at the right level of detail
for first use:

1. Start from a cleaned, unmodified protein PDB.
2. Define each attachment as a site, moiety, and mechanism.
3. Use `build_conjugate()` for direct request objects or
   `build_conjugate_from_config()` for config-driven NHS-Lys polymer builds.
4. Inspect `ConjugationResult` paths before moving on to simulation.
5. Use built-in mechanisms first, and treat custom SMARTS as a validation and
   extension design tool until a concrete mechanism implementation exists.

For working examples, see the notebooks in
`src/polyzymd/builders/conjugation/poc/`, especially
`public_api_conjugation_walkthrough.ipynb` for NHS-Lys and
`public_api_n_glycosylation_walkthrough.ipynb` for N-glycosylation.
