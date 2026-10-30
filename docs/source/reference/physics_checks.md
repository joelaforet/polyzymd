# Physics checks

Every PolyzyMD build runs each of these checks on the system the requested engine
will run (OpenMM `system.xml`, or GROMACS `.top`/`.gro` with the production `.mdp`),
and records the result in `build_provenance.json` and `build_provenance.md`.
A check whose default is FAIL stops the build; WARN and INFO are recorded only.

To keep a deliberate departure, accept the failing (or warning) check in the config with a
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

| Check | Runs on | Default | Catches |
|---|---|---|---|
| `charges.site_integer` | each component (protein) | FAIL | a modification site whose charges, with the rest of its residues, do not sum to their formal charge |
| `charges.site_shift` | each component (protein) | FAIL | a site whose NAGL charges need a uniform shift above the bound to reach its integer charge (the capped model and the site disagree) |
| `charges.source` | each component | FAIL | partial charges that differ from the declared source (library, NAGL model, preset file) |
| `components.declared` | the system | FAIL | a component, its segment ID, force field or charge source that differs from the config |
| `composition.counts` | the system | FAIL | a molecule or chain count that differs from the config |
| `composition.density` | the system | INFO | the density the box reached, against the config's target |
| `composition.polymer_mass_fraction` | the system (polymer) | INFO | the polymer weight fraction of polymer plus solvent the build reached |
| `composition.polymer_sequences` | the system (polymer) | INFO | the polymer sequences used, and the monomer fractions they reach |
| `constraints.hbonds` | the system | FAIL | a bond to hydrogen left flexible at a time step that needs it constrained, or a step above 2 fs without hydrogen mass repartitioning |
| `engine.realized` | the system | FAIL | GROMACS run parameters (.mdp) other than the config's: integrator, thermostat, time step, length, temperature or annealing, pressure coupling, bond constraints, velocity generation or position-restraint switches |
| `geometry.buried_water` | the system (protein) | WARN | water buried inside the protein (may be deliberate, as a crystal water) |
| `geometry.image_bridging` | the system (polymer) | FAIL | a polymer chain touching both the protein and the protein's periodic image |
| `geometry.interchain_overlap` | the system | WARN | polymer chains clashing with each other (usually resolved by minimization; a sign of a forced Packmol result) |
| `geometry.modified_protein_clash` | each component (protein) | FAIL | atoms of a modified protein without a bond or angle between them closer than 1.0 Å (heavy-heavy, hydrogen-heavy) or 0.8 Å (hydrogen-hydrogen): a fragment or side chain placed through the protein, which minimization can only fix by stretching bonds or moving the protein |
| `geometry.periodic_image` | the system | FAIL | a fixed component closer to its own periodic image than the box padding guarantees |
| `geometry.polymer_overlap` | the system | FAIL | a polymer heavy atom clashing with the protein or ligand |
| `geometry.water_overlap` | the system | FAIL | a solvent heavy atom clashing with a solute heavy atom (closer than half the Packmol tolerance the build used) |
| `ions.ionic_strength` | the system | FAIL | an ionic strength more than one ion away from the configured target |
| `ions.neutral` | the system | FAIL | a net charge left on a system the config asked to neutralize |
| `ions.salt` | the system | FAIL | NaCl pairs more than one pair away from the configured concentration |
| `masses` | each component | FAIL | atom masses that differ from the force field's (unintended hydrogen mass repartitioning) |
| `modeling_exceptions.unused` | the system | WARN | a modeling_exceptions entry that accepted nothing (stale in a published config) |
| `modifications.applied` | each component (protein) | FAIL | a planned modification missing from the built system: a bond it forms, an atom it adds or removes, or an atom renamed |
| `modifications.chemistry` | each component (protein) | FAIL | a modified protein whose chemistry (elements, bond orders, formal charges) differs from an independent reading of its structure with openff-pablo |
| `modifications.declared` | each component (protein) | FAIL | a built protein whose modifications or parameter sets differ from the config's mutations:, modifications: and parameter sets |
| `modifications.protonation_state` | each component (protein) | WARN | a mutation that leaves a side chain in a state other than the usual one at pH 7 (a CCD component's neutral Asp or Glu, or doubly protonated His) |
| `modifications.sites` | each component (protein) | FAIL | a modification site whose atoms are not those its changes define (side chain from CB of the anchor, touched backbone atoms, added atoms, attached residues), or that differ from the component record's |
| `nonbonded.settings` | the system | FAIL | a cutoff, switching distance, PME or dispersion correction other than the force field's |
| `parameters.bonded` | each component | FAIL | bonds, angles, torsions or impropers not from the component's declared force field |
| `parameters.constraints` | each component | FAIL | constraints missing, extra, or at lengths other than the declared force field's |
| `parameters.coverage` | the system | FAIL | atoms no component owns, molecules other than the components declare, or terms (bonds, constraints, exceptions) joining two molecules |
| `parameters.glycan_force_field` | each component (protein) | WARN | a glycan not parameterized with GLYCAM06j, the recommended glycan force field (general force fields give stiffer glycosidic angles and other torsion minima) |
| `parameters.modified_backbone` | each component (protein) | FAIL | the backbone (N, H, CA, HA, C, O) of a modified residue not on ff14SB: its charges and terms must match its unmodified neighbours' |
| `parameters.nonbonded` | each component | FAIL | Lennard-Jones parameters not from the component's declared force field |
| `parameters.one_four` | each component | FAIL | 1-2 and 1-3 exclusions or 1-4 pairs not as the owning force field scales them (GLYCAM06j leaves some glycan 1-4 pairs unscaled; one global scale gets them wrong) |
| `parameters.rosemary_alpha` | each component (protein) | WARN | a protein charged whole with NAGL by its force field (OpenFF Rosemary): an alpha release is not for production, and OpenFF pairs Rosemary with OPC3 water; reported only for such proteins |
| `restraint.atom_identity` | the system | FAIL | a restrained atom that is not the atom its selection names in its component |
| `restraint.parameters` | the system | FAIL | a restraint whose form, distance or force constant differs from the config, or a restraint missing or extra |
