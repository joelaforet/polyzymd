# Analyze the trajectories with your own tools

PolyzyMD writes standard files: an Amber topology and DCD trajectories for
OpenMM, and a GROMACS topology, run input and XTC trajectories for GROMACS.
MDAnalysis, MDTraj, VMD and the `gmx` tools read them directly. This page
says which files to read, which names the built system uses, how the frame
times are laid out, and how to bring your results back into a project.

The examples use a project with three conditions of Trp-cage: `Water` and
`Urea 1 M`, run with OpenMM, and `GROMACS Water`, run with GROMACS. The Python
examples run from the project folder.

## Find the replicate folder

Each replicate runs in its own folder: the config's `naming_template` (such as
`trpcage_300K_run1`) under its `scratch_directory`, or under its
`projects_directory` when it has no scratch directory. In the example project
the folders are:

```text
runs/trpcage/water/trpcage_300K_run1/            # OpenMM
runs/trpcage/urea_1_m/trpcage_300K_run1/         # OpenMM
runs/trpcage/gromacs_water/trpcage_300K_run1/    # GROMACS, with a gromacs/ folder
```

## Which files to read

### OpenMM replicate

| File | Holds | Use it for |
|---|---|---|
| `system.prmtop` | Every atom, residue, element, mass, charge and bond. One segment `SYSTEM`, no chain IDs | The topology of your analysis |
| `production_N/production_N_trajectory.dcd` | The production frames of segment `N` | The coordinates. Join the segments in the order of `N` |
| `solvated_system.pdb` | The built coordinates, with segment IDs (`PROT`, `SOL`, `ION`, ...) and chain IDs (`A`, `D`, ...). No bonds and no charges | Chain IDs and segment IDs for your topology, and viewers such as PyMOL and VMD |
| `production_N/production_N_state_data.csv` | Step, time in ps, energies, temperature, volume and density of each frame | Frame times and thermodynamic quantities |
| `progress.json` | One record per equilibration stage and production segment | Which segments are complete, and how they join |

Do not use `production_N/production_N_topology.pdb` or the files of
`equilibration_*/` as the topology of an analysis. `production_N_topology.pdb`
has the same atoms in the same order, but OpenMM wrote its chain IDs: they
run from `A` to `Z` in blocks of a few hundred atoms and do not mark the
protein or the solvent.

### GROMACS replicate

The run files are in `gromacs/`, named after the config's prefix (here
`trpcage`). `solvated_system.pdb` is in the replicate folder.

| File | Holds | Use it for |
|---|---|---|
| `gromacs/prod.tpr` | The compiled production input: atoms, bonds, masses, charges | The topology for `gmx` tools, and for MDAnalysis when it can read the file version |
| `gromacs/trpcage.top` and its `.itp` files | The same, as text | The topology for MDAnalysis when it cannot read `prod.tpr` |
| `gromacs/prod_centered.xtc` | Production frames with whole molecules and the protein centred | Most analyses. PolyzyMD's own analyses read it first |
| `gromacs/prod_nojump.xtc` | Production frames where no atom jumps across the box between frames | Displacements over time. Molecules split by the box in the first frame stay split |
| `gromacs/prod.xtc` | Production frames as `mdrun` wrote them, each atom in the box | The raw record. Molecules can be split across the box |
| `gromacs/components.ndx` | One index group per component, under its name and its segment ID | Selections for `gmx` tools |
| `gromacs/progress.json` | One record per stage | Which stages are complete |
| `solvated_system.pdb` | Segment IDs and chain IDs, as for OpenMM | Chain IDs and segment IDs for your topology |

MDAnalysis 2.10 reads TPR files of GROMACS 2025 and older, not those of
GROMACS 2026. The example below reads the `.top` file instead.

## Names in the built system

The same names are used in `system.prmtop`, `solvated_system.pdb` and the
GROMACS topology.

| Component | Residue names | Atom names | Segment ID | Chain ID |
|---|---|---|---|---|
| Protein | Standard residue names, numbered from 1 | PDB names (`N`, `CA`, `C`, `O`, `CB`, ...) | `PROT` | `A` |
| Small molecule (ligand) | `residue_name` of the solute, `LIG` by default | The SDF's atom names; without them, element, count and `x`: `C1x`, `C2x`, `O1x`, `H1x` | `LIG` | `B` |
| Polymer | `residue_name` of each monomer, by default the first three letters of its `name`, upper case | Set by the polymer's SDF; list them as shown below | `POLY` | `C` |
| Water | `HOH` | `O1x`, `H1x`, `H2x` | `SOL` | `D`, then `E`, ... |
| Ions | `Na+`, `Cl-` | `NA`, `CL` | `ION` | `D`, then `E`, ... |
| Co-solvents | Their residue name (urea `URE`, DMSO `DMS`; see the co-solvent library in {doc}`../reference/configuration`) | Element, count and `x`: urea has `C1x`, `O1x`, `N1x`, `N2x`, `H1x` to `H4x` | `SOL` | `D`, then `E`, ... |

A solute's segment ID is the `segid` you set in its config, if you set one.
Water, ions and co-solvents fill chain `D`. Their residues are numbered from 1
and continue on chain `E` after residue 9999. In `system.prmtop` the same
residues are numbered on from the solutes instead. Select solvent by residue
name or segment ID, not by residue number.

Usual selections such as `name OW`, `resname SOL`, `resname NA CL` or
`resname HOH and name O` match no atoms in these systems. These MDAnalysis
selections work:

| Atoms | Selection | Needs the IDs of `solvated_system.pdb` |
|---|---|---|
| Protein | `protein`, or `segid PROT`, or `chainID A` | Only for `segid` and `chainID` |
| Water oxygens | `resname HOH and name O1x` | No |
| Ions | `resname Na+ Cl-`, or `segid ION` | Only for `segid` |
| Water and co-solvents | `segid SOL` | Yes |
| Urea | `resname URE` | No |
| Heavy atoms of the protein | `protein and not element H` | For GROMACS, which has no elements |

List the residue and atom names of any built system with MDAnalysis:

```python
import MDAnalysis as mda

u = mda.Universe("runs/trpcage/urea_1_m/trpcage_300K_run1/system.prmtop")
for resname in sorted(set(u.residues.resnames)):
    residue = u.select_atoms(f"resname {resname}").residues[0]
    print(resname, " ".join(residue.atoms.names))
```

## Time and frames

- The trajectories store time in ps. MDAnalysis reports times in ps and
  lengths in Å. PolyzyMD's configs, `study.yaml` windows such as
  `equilibration: 10ns`, and the `times=` argument of your own functions use
  ns.
- Production time starts at 0 at the start of production. Equilibration is not
  in the production trajectories.
- OpenMM writes its first frame one report interval after the start. A
  production of `duration: 0.03` (ns) with `samples: 30` gives a DCD of 30
  frames, at 1, 2, ... 30 ps.
- GROMACS also writes the frame at 0 ps. The same settings give 31 frames, at
  0, 1, ... 30 ps.
- The study's equilibration window is not applied when you read the files
  yourself. Leave out the frames before it.

### Join OpenMM segments

A run that was continued or resumed has several segments: `production_0/`,
`production_1/`, and so on. Sort them by the number `N`, not by name:
`production_10` sorts before `production_2` as text. MDAnalysis joins a list
of DCD files into one trajectory:

```python
import re
from pathlib import Path

import MDAnalysis as mda

rep = Path("runs/trpcage/urea_1_m/trpcage_300K_run1")


def segment_number(path):
    return int(re.search(r"production_(\d+)_trajectory", path.name).group(1))


segments = sorted(rep.glob("production_*/production_*_trajectory.dcd"), key=segment_number)
u = mda.Universe(rep / "system.prmtop", segments)
print(u.trajectory.n_frames, u.trajectory[0].time, u.trajectory[-1].time)
```

The segments of a chain written by PolyzyMD 1.3.0 follow on from one another:
each starts one report interval after the previous one ends, and the DCD
headers carry the real times. A segment folder that holds an `INTERRUPTED`
file belongs to the chain too; its frames are complete.

`progress.json` records how each segment joins the one before it:

- `overlap_frames`: frames at the start of the segment that the previous
  segment had already written. PolyzyMD 1.3.0 does not write such frames
  again, so this is 0.
- `gap_frames`: report steps between the two segments that no segment wrote.
  The frames after a gap keep their times.

Segments written by older PolyzyMD versions can have two faults:

- A segment starts at the time of the previous segment's last frame. Keep one
  frame per time. PolyzyMD's loader keeps the later one:

  ```python
  import numpy as np

  times = np.array([ts.time for ts in u.trajectory])
  step = np.rint(times / u.trajectory.dt).astype(int)
  keep = np.flatnonzero(np.r_[step[1:] != step[:-1], True])
  for ts in u.trajectory[keep]:
      pass  # your analysis
  ```

- Under OpenMM 8.2 or later, the DCD header of a continued segment gives the
  wrong first step, and its times start again from one report interval. Take
  the times from the state-data files instead:

  ```python
  import pandas as pd

  csv_times = np.concatenate([
      pd.read_csv(p.with_name(p.name.replace("_trajectory.dcd", "_state_data.csv")))["Time (ps)"]
      for p in segments
  ])
  ```

A GROMACS run has no segments: a continued run appends to `prod.xtc`, and
`gmx trjconv` writes `prod_nojump.xtc` and `prod_centered.xtc` from it.

## Load a replicate in MDAnalysis

### Add chain IDs and segment IDs

`system.prmtop` has one segment, `SYSTEM`, and no chain IDs. The GROMACS
topology names its segments after molecule types (`MOL0`, `water_TIP3P`,
`MOL1942`, ...) and has no elements. `solvated_system.pdb` has the atoms in
the same order. Copy its chain IDs, segment IDs and, where missing, elements
onto your universe:

```python
def add_build_ids(u, pdb_file):
    """Copy chain IDs, segment IDs and elements from solvated_system.pdb."""
    pdb = mda.Universe(pdb_file)
    assert (pdb.atoms.names == u.atoms.names).all(), "atoms differ from the build"
    u.add_TopologyAttr("chainIDs", pdb.atoms.chainIDs)
    if not hasattr(u.atoms, "elements"):
        u.add_TopologyAttr("elements", pdb.atoms.elements)
    segids = pdb.residues.segids
    for segid in dict.fromkeys(segids):
        u.residues[segids == segid].segments = u.add_Segment(segid=segid)
```

You can also load `solvated_system.pdb` as the topology. It has segment IDs,
chain IDs and elements, but no bonds and no charges, and its solvent residue
numbers differ from those of `system.prmtop`.

### OpenMM: radius of gyration and contacts

With the universe `u` of the segments above:

```python
import numpy as np
from MDAnalysis.lib.distances import capped_distance

add_build_ids(u, rep / "solvated_system.pdb")
protein = u.select_atoms("segid PROT")
heavy = u.select_atoms("protein and not element H")
urea = u.select_atoms("resname URE and not element H")

rg, n_urea = [], []
for ts in u.trajectory:
    rg.append(protein.radius_of_gyration())
    pairs = capped_distance(heavy.positions, urea.positions, max_cutoff=4.0,
                            box=ts.dimensions, return_distances=False)
    n_urea.append(len(np.unique(urea[pairs[:, 1]].resindices)))
print(f"Rg {np.mean(rg):.2f} A, urea molecules within 4 A: {np.mean(n_urea):.2f}")
```

### GROMACS: radius of gyration and contacts

Continuing from the code above, read the `.top` file with
`topology_format="ITP"` and `infer_system=True`, so that MDAnalysis builds
every molecule that `[ molecules ]` lists:

```python
gmx = Path("runs/trpcage/gromacs_water/trpcage_300K_run1/gromacs")
g = mda.Universe(gmx / "trpcage.top", gmx / "prod_centered.xtc",
                 topology_format="ITP", infer_system=True)
add_build_ids(g, gmx.parent / "solvated_system.pdb")

protein = g.select_atoms("segid PROT")
heavy = g.select_atoms("protein and not element H")
ions = g.select_atoms("segid ION")

rg, n_ions = [], []
for ts in g.trajectory:
    rg.append(protein.radius_of_gyration())
    pairs = capped_distance(heavy.positions, ions.positions, max_cutoff=4.0,
                            box=ts.dimensions, return_distances=False)
    n_ions.append(len(np.unique(pairs[:, 1])))
print(f"Rg {np.mean(rg):.2f} A, ions within 4 A: {np.mean(n_ions):.3f}")
```

### Periodic boundaries

- The box can be triclinic, such as a rhombic dodecahedron. Pass
  `box=ts.dimensions` to distance functions.
- OpenMM DCD frames keep each molecule whole, wrapped into the box as a
  unit. A molecule can jump across the box between frames. The radius of
  gyration and contacts computed with `box=` are not affected. Displacements
  over time, such as a mean squared displacement, need the jumps removed,
  for example with MDAnalysis's `NoJump` transformation.
- In GROMACS, use `prod_centered.xtc` unless you need displacements over time.

## Use `gmx` tools

Run the tools in the `gromacs/` folder of the replicate. `prod.tpr` is the
structure file. `components.ndx` has the groups `trpcage` and `PROT` (the
protein, under its name and its segment ID), `solvent` and `SOL` (water and
co-solvents), and `ions` and `ION`.

```bash
cd runs/trpcage/gromacs_water/trpcage_300K_run1/gromacs
gmx gyrate -s prod.tpr -f prod_centered.xtc -n components.ndx -sel PROT -o gyrate.xvg
gmx select -s prod.tpr -f prod_centered.xtc -n components.ndx \
    -select 'group ION and within 0.4 of (group PROT and not name "H.*")' -os ion_count.xvg
```

`gmx` tools report lengths in nm. In a `gmx select` selection, a quoted
string with characters such as `+` is a regular expression, so
`resname "Na+"` matches no atoms. Select ions by `name NA CL` or
`group ION`. `gmx make_ndx` takes the residue names as they are:

```bash
printf 'r HOH & a O1x\nr Na+ | r Cl-\nq\n' | gmx make_ndx -f prod.tpr -n components.ndx -o mine.ndx
```

## Bring your results into the project

### As a `function:` analysis

Write the measurement as a Python function of AtomGroups and list it under
`analyses:` in `project.yaml` (or `study.yaml`). PolyzyMD then runs it on every
replicate of every study, over the study's equilibration window, and stores,
compares and reports it like a shipped analysis. Its loader joins the segments,
adds the chain IDs, segment IDs and elements, and reads GROMACS 2026 runs.
`polyzymd study freeze` and `polyzymd project freeze` warn when the function's
code changed after its results were stored.

```python
# analyses/mine.py
import numpy as np
from MDAnalysis.lib.distances import capped_distance


def my_rg(protein):
    """Mass-weighted radius of gyration of the protein, in A."""
    return float(protein.radius_of_gyration())


def cosolvent_contacts(protein, cosolvent):
    """Co-solvent molecules with a heavy atom within 4 A of a protein heavy atom."""
    if len(cosolvent) == 0:
        return 0.0
    pairs = capped_distance(protein.positions, cosolvent.positions, 4.0,
                            box=protein.dimensions, return_distances=False)
    return float(len(np.unique(cosolvent[pairs[:, 1]].resindices)))
```

```yaml
analyses:
  my_rg:
    function: analyses/mine.py:my_rg
    kind: timeseries
    unit: A
    selections: {protein: protein}
  my_urea_contacts:
    function: analyses/mine.py:cosolvent_contacts
    kind: timeseries
    unit: count
    allow_empty: true          # the Water condition has no urea
    selections: {protein: protein and not element H, cosolvent: resname URE and not element H}
```

```bash
polyzymd project check .
polyzymd analyze --project .
```

{doc}`study_yaml` lists every key of a `function:` entry.

### As a script in `stats/`

Results that a `function:` entry cannot produce, such as the output of `gmx`
tools, go in `stats/` with the script that made them. `polyzymd project
freeze` then treats them as follows:

- It deposits every file under `stats/`, and hashes each one (size and
  SHA-256) in the manifest.
- It does not record which trajectories or which version of the script made
  an output file, and it does not check whether the output is older than the
  script. Run the script again before you freeze.
- Freeze stops while a file in `stats/` is not committed. Commit the script
  and its outputs.
- The deposit's `README.md` says how to redraw the figures and recompute
  PolyzyMD's analyses. It does not mention your scripts. Say how to run them
  in `stats/README.md`.
- Files in folders that freeze does not know, such as `own_analysis/`, are
  not deposited; freeze names them in a warning.

Do not put files of your own in `results/`. Freeze deposits and hashes them
with PolyzyMD's results, with no warning, and `polyzymd analyze --recompute`
does not reproduce them.

## Related pages

- {doc}`study_yaml` — run your own function in a study
- {doc}`study_api` — get a replicate's `Universe` from PolyzyMD's loader
- {doc}`project` — statistics scripts and freezing a project
- {doc}`../reference/data_requirements` — the files of a replicate folder
- {doc}`../explanation/residue_assignment` — chain IDs and residue numbers
