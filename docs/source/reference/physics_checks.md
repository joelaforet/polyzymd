# Physics checks

Every PolyzyMD build runs each of these checks on the system the requested engine
will run (OpenMM `system.xml`, or GROMACS `.top`/`.gro` with the production `.mdp`),
and records the result in `build_provenance.json` and `build_provenance.md`.
A check whose default is FAIL stops the build; WARN and INFO are recorded only.

To keep a deliberate departure, accept the failing check in the config with a
reason. The build record then lists the accepted finding, its reason and the
atoms it covered:

```yaml
modeling_exceptions:
  - check: geometry.buried_water
    component: protein        # optional; every component the check runs on if left out
    reason: "Crystal waters W301 and W305 kept deliberately (PDB 4K6G)"
```

An exception that matches no failing check is reported, so stale entries do not
linger in a published config.

This page is generated from PolyzyMD's check registry
(`polyzymd.audit.registry.catalog_markdown`).

| Check | Runs on | Default | Catches | Origin |
|---|---|---|---|---|
| `charges.source` | each component | FAIL | partial charges that differ from the declared source (library, NAGL model, preset file) | LipA production runs that recomputed charges |
| `components.declared` | the system | FAIL | a component, its segment ID, force field or charge source that differs from the config | parameter checks compare with the builder's records; this ties the records to the config |
| `composition.counts` | the system | FAIL | a molecule or chain count that differs from the config |  |
| `composition.density` | the system | INFO | the density the box reached, against the config's target |  |
| `composition.polymer_mass_fraction` | the system (polymer) | INFO | the polymer weight fraction of polymer plus solvent the build reached | 20 wt% decision for the Paper 1 rerun |
| `composition.polymer_sequences` | the system (polymer) | INFO | the polymer sequences used, and the monomer fractions they reach |  |
| `constraints.hbonds` | the system | FAIL | a bond to hydrogen left flexible at a time step that needs it constrained, or a step above 2 fs without hydrogen mass repartitioning |  |
| `geometry.buried_water` | the system (protein) | WARN | water buried inside the protein (may be deliberate, as a crystal water) | legacy CALB/RML builds |
| `geometry.image_bridging` | the system (polymer) | FAIL | a polymer chain touching both the protein and the protein's periodic image |  |
| `geometry.interchain_overlap` | the system | WARN | polymer chains clashing with each other (usually resolved by minimization; a sign of a forced Packmol result) | LipA July builds |
| `geometry.periodic_image` | the system | FAIL | a fixed component closer to its own periodic image than the box padding guarantees | old brick boxes |
| `geometry.polymer_overlap` | the system | FAIL | a polymer heavy atom clashing with the protein or ligand |  |
| `geometry.water_overlap` | the system | FAIL | a solvent heavy atom clashing with a solute heavy atom (closer than half the Packmol tolerance the build used) | legacy CALB/RML builds with water inside the protein |
| `ions.ionic_strength` | the system | FAIL | an ionic strength more than one ion away from the configured target | nacl_concentration meaning change; buffer matched by ionic strength (Paper 1 rerun) |
| `ions.neutral` | the system | FAIL | a net charge left on a system the config asked to neutralize |  |
| `ions.salt` | the system | FAIL | NaCl pairs more than one pair away from the configured concentration | nacl_concentration meaning change |
| `masses` | each component | FAIL | atom masses that differ from the force field's (unintended hydrogen mass repartitioning) |  |
| `nonbonded.settings` | the system | FAIL | a cutoff, switching distance, PME or dispersion correction other than the force field's |  |
| `parameters.bonded` | each component | FAIL | bonds, angles, torsions or impropers not from the component's declared force field | 1.3 force-field order bug: Sage terms on ff14SB proteins (Paper 1) |
| `parameters.constraints` | each component | FAIL | constraints missing, extra, or at lengths other than the declared force field's | 1.3 force-field order bug: Sage terms on ff14SB proteins (Paper 1) (Sage bond-to-hydrogen lengths on the protein) |
| `parameters.coverage` | the system | FAIL | atoms no component owns, molecules other than the components declare, or terms (bonds, constraints, exceptions) joining two molecules | every term must have exactly one declared source |
| `parameters.nonbonded` | each component | FAIL | Lennard-Jones or 1-4 terms not from the component's declared force field | 1.3 force-field order bug: Sage terms on ff14SB proteins (Paper 1) |
| `restraint.atom_identity` | the system | FAIL | a restrained atom that is not the atom its selection names in its component | pdbindex off-by-one in Paper 1 restraint selections |
| `restraint.parameters` | the system | FAIL | a restraint whose form, distance or force constant differs from the config, or a restraint missing or extra |  |
