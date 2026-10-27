"""Physics audit: proof that a built system is the physics its config asked for.

Every build is checked against the components that made it
(:mod:`polyzymd.core.components`): each component is rebuilt on its own from the
force field and charge source the config names, and every term of the system
the requested engine will run is compared with that reference. System-wide and
geometric checks follow. This package never imports the builders, so it cannot
share a mistake with the code it checks.

- :mod:`~polyzymd.audit.registry`: the catalog of checks.
- :mod:`~polyzymd.audit.realized`: the system an engine will run, read from its files.
- :mod:`~polyzymd.audit.terms`: every parameter of a System, indexed for comparison.
- :mod:`~polyzymd.audit.molecules`: matching atoms to reference molecules, and
  single-force-field reference builds.
"""
