# Run a PolyzyMD build with your own protocol

PolyzyMD builds a system, audits it, and runs a standard protocol: minimization,
equilibration stages and production. Some simulations need more: enhanced
sampling with PLUMED, your own GROMACS `.mdp` files, or an OpenMM script of
your own. Run those yourself, then bring the trajectories back into the study,
where `polyzymd analyze` reads them like any other replicate.

PolyzyMD then vouches for the system it built, which the physics audit checked
and `build_provenance.json` records. It does not vouch for the run: the
protocol is yours, and the replicate's records say so.

The steps are the same for both engines:

1. Build the replicate with PolyzyMD.
2. Put the inputs of your run in `protocol/`.
3. Run it with your protocol, in a folder of its own.
4. Copy the trajectory back to where analysis reads it.
5. Analyze.

## 1. Build

Build the replicate as usual. The build runs the physics audit; read
`build_provenance.md` before you go on.

```bash
polyzymd build -c conditions/water/config.yaml -r 1                     # OpenMM
polyzymd build -c conditions/water/config.yaml -r 1 --format gromacs    # GROMACS
```

Do not run `polyzymd run` or `polyzymd submit` on this replicate.

## 2. Declare the protocol, before you run

Make the folder `protocol/` in the replicate folder for OpenMM, or
`gromacs/protocol/` for GROMACS, and put the inputs of your run there: your
scripts, `.mdp` files and PLUMED input. Make it before you run, so PolyzyMD
leaves the replicate alone while your run is going. After the run, add any
output a reader needs to judge it, such as `COLVAR` and `HILLS`.

```bash
mkdir protocol && cp run_metad.py protocol/                 # OpenMM
mkdir gromacs/protocol && cp minim.mdp nvt.mdp md.mdp gromacs/protocol/   # GROMACS
```

The folder is how PolyzyMD knows the replicate is run outside it. With it:

- `polyzymd run`, `submit`, `recover` and `build` refuse the replicate, and a
  PolyzyMD job already queued for it stops without running, so your files are
  never overwritten or deleted;
- `polyzymd status` shows the replicate as run with its own protocol, not as
  dead, and suggests no resubmission;
- each analysis record of the replicate holds `run_protocol: own`;
- a frozen study deposits the files of `protocol/` with the built topology
  (and `prod.tpr` for GROMACS), never the `.mdp` files PolyzyMD wrote but did
  not run; its manifest records the run as `own protocol`, and its checklist
  marks the replicate's settings as its own.

Without the folder, a frozen study warns that it cannot say how the
replicate was run.

A GROMACS job script or `run_<prefix>_gromacs.sh` written by a PolyzyMD older
than this check runs PolyzyMD's own protocol whatever the folder holds: do
not run one in this replicate, and cancel any of its jobs still queued.

## 3. Run it with your protocol

Run in a new folder inside the replicate folder, such as `own/`, so the built
files stay as they were audited.

**OpenMM.** Load the built `system.xml` and `solvated_system.pdb`, and add
what your protocol needs. This script runs metadynamics on the radius of
gyration of the protein's Cα atoms with the `openmm-plumed` plugin:

```python
import openmm as mm
import openmm.app as app
from openmm import unit
from openmmplumed import PlumedForce

pdb = app.PDBFile("../solvated_system.pdb")
with open("../system.xml") as f:
    system = mm.XmlSerializer.deserialize(f.read())
ca = [a.index + 1 for a in pdb.topology.atoms()
      if a.name == "CA" and a.residue.chain.id == "A"]   # chain A is the protein
system.addForce(PlumedForce(f"""
rg: GYRATION ATOMS={",".join(map(str, ca))}
metad: METAD ARG=rg SIGMA=0.02 HEIGHT=1.0 PACE=100 FILE=HILLS
PRINT ARG=rg,metad.bias STRIDE=100 FILE=COLVAR
"""))
integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                         2 * unit.femtoseconds)
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
simulation.reporters.append(app.DCDReporter("traj.dcd", 100))
simulation.step(10000)
```

Save it as `protocol/run_metad.py` and run it from `own/`:

```bash
mkdir own && cd own
python ../protocol/run_metad.py
cp COLVAR HILLS ../protocol/
```

**GROMACS.** Use the exported topology and coordinates in `gromacs/` with
your own `.mdp` files:

```bash
cd gromacs
mkdir own
gmx grompp -f protocol/minim.mdp -c trpcage.gro -p trpcage.top -o own/em.tpr
gmx mdrun -deffnm own/em
gmx grompp -f protocol/nvt.mdp -c own/em.gro -p trpcage.top -o own/nvt.tpr
gmx mdrun -deffnm own/nvt
gmx grompp -f protocol/md.mdp -c own/nvt.gro -t own/nvt.cpt -p trpcage.top -o own/md.tpr
gmx mdrun -deffnm own/md
```

Use the exported `.top` and `.itp` files as they are. If your protocol
changes the system itself, such as a bias, an extra restraint or a changed
parameter, that change is part of your protocol: keep its input in
`protocol/`.

## 4. Copy the trajectory back

Analysis finds trajectories by name.

**OpenMM.** Copy the trajectory to `production_0/production_0_trajectory.dcd`
in the replicate folder. A run in several parts goes to `production_1/`,
`production_2/` and so on, with no gaps; their times must follow on from one
another. Analysis reads the topology from the built `system.prmtop`.

```bash
mkdir production_0
cp own/traj.dcd production_0/production_0_trajectory.dcd
```

**GROMACS.** Copy the production run input, trajectory and log into
`gromacs/` as `prod.tpr`, `prod.xtc` and `prod.log`. Then make the molecules
whole and centre the protein, as PolyzyMD does after its own runs; otherwise
a protein split across the periodic boundary gives wrong radii of gyration,
distances and contacts. Analysis reads `prod_centered.xtc` first, and warns
when it can only read the raw `prod.xtc`.

```bash
cd gromacs
cp own/md.tpr prod.tpr
cp own/md.xtc prod.xtc
cp own/md.log prod.log      # the deposit reads the GROMACS version from it
echo "System" | gmx trjconv -s prod.tpr -f prod.xtc -o prod_nojump.xtc -pbc nojump
echo -e "Protein\nSystem" | gmx trjconv -s prod.tpr -f prod_nojump.xtc \
    -o prod_centered.xtc -center -pbc mol -ur compact
```

## 5. Analyze

```bash
polyzymd status -c conditions/water/config.yaml --no-slurm
polyzymd analyze --project .
```

The replicate takes part in every analysis and comparison like the
replicates PolyzyMD ran. Frame times come from the trajectory, and the
study's equilibration window applies to them. To read the trajectories with
MDAnalysis or `gmx` tools yourself, see {doc}`own_analysis`.

```{warning}
PolyzyMD analyses weight every frame equally. The frames of a biased run, such
as metadynamics or umbrella sampling, do not sample the unbiased ensemble:
reweight them before you report an average from them. For metadynamics the
bias changes during the run, so the weights need its time-dependent
correction (in PLUMED, `CALC_RCT` and the `rbias` it gives), not the bias
alone. PolyzyMD does not reweight.
```

## Run the replicate with PolyzyMD instead

Move `protocol/` and the trajectories (`production_N/`, or the `prod*` files
in `gromacs/`) out of the replicate folder. It is then a built replicate that
has not started, and `polyzymd run` or `submit` runs it with PolyzyMD's
protocol.

## Start from the end of a finished replicate

PolyzyMD cannot yet start a new condition from the final state of a finished
replicate, for example to heat it after 1 µs at 300 K. To do it anyway, build
the new condition's replicate as above, put your protocol in its `protocol/`
folder, and start that protocol from the finished replicate's last
`production_N_state.xml`. Analysis then reads it as a run with your own
protocol. To keep PolyzyMD's own runs, build and equilibrate the new
condition from its own config instead.

## Related pages

- {doc}`../explanation/physics_audit` — what the build audit checks
- {doc}`run_gromacs` — run GROMACS with PolyzyMD's protocol
- {doc}`study_folder` — studies, conditions and replicates
- {doc}`own_analysis` — analyze the trajectories with your own tools
