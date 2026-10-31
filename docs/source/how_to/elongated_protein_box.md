# Choose a box for an elongated protein

```{index} !rectangular box, !per-axis padding, !box shape
```

Size the periodic box for a long, thin protein (a multi-domain fibronectin
fragment, an antibody, a fibre) without filling most of it with water.

## Why the default box is large

A cube or a rhombic dodecahedron has one edge: the solute's largest
atom-to-atom distance plus `2 * padding`. That keeps the solute
`2 * padding` from its periodic images in any orientation, so the protein may
turn freely during the run. For a protein three times longer than it is wide,
almost all of that box is water. A 4-domain fibronectin fragment in a rhombic
dodecahedron with 1 nm padding is about 640,000 atoms.

## Step 1: Use a rectangular box padded per axis

```yaml
solvent:
  box:
    shape: rectangular
    padding: {x: 1.2, y: 1.2, z: 1.5}   # nm from the protein to the box edge, per axis
```

Before sizing the box, PolyzyMD turns the fixed solutes (the protein with
everything attached to it, and any ligand) together so that their longest
principal axis lies on z and the shortest on x. Each box edge is then the
protein's extent along that axis plus twice that axis's padding. A single
number (`padding: 1.2`) pads every axis the same.

The rotation comes from the mass-weighted heavy atoms of the protein segment
and any ligand. A modified protein's attachments (`G###`, `P###`, `M###`
segments) are left out, because grafted chains can differ between
replicates and every replicate of a condition must share one box. When two
principal moments are within 1%, the axes in their plane are not defined by
the shape, so that plane keeps the input orientation (the build log says so).
The rotation is recorded in `build_provenance.json` (`box_alignment`), and
`modified_protein.pdb` is written in the rotated frame, matching the solvated
system. To keep the input orientation, add `align: none`.

## Step 2: Respect the minimum padding

Every axis needs at least the nonbonded cutoff plus 0.1 nm (1.0 nm with the
0.9 nm cutoff of Sage and the ff14SB port). `validate` and `build` refuse a
smaller padding and name the axis. Pad the long axis more if the protein is
flexible along it.

## Step 3: Check the protein's distance to its own image

A rectangular box is sized for the starting orientation. A long protein that
turns during the run can come within the cutoff of its own periodic image,
which then interacts with it. The build audit reminds you of this with a
warning (`geometry.rectangular_box`). After the run, measure the closest
approach:

```python
import MDAnalysis as mda

from polyzymd.analyses.periodic_images import image_distance_findings, solute_image_distances

u = mda.Universe("solvated_system.pdb", "production.dcd")
distances = solute_image_distances(u, "segid PROT and not name H*", step=10)
for warning in image_distance_findings(distances, nonbonded_cutoff_nm=0.9):
    print(warning)
```

The protein is treated as one piece: its first frame must be whole (it is in
a PolyzyMD build), and later frames are unwrapped with MDAnalysis' `NoJump`,
so chains in contact across the box edge are not counted as contacts with
the image. No bonds or masses are needed.

If the protein comes within the cutoff of its image, rebuild with more
padding, or use a rhombic dodecahedron, which is sized for any orientation.
