# Replicate records

`MDAUniversePolicy` records which files one replicate was measured from: the
condition label, the replicate number, and the provenance the trajectory loader
already resolved. Every replicate artifact carries it, so a stored number can
name the topology and the trajectory segments behind it.

`MDAJobResult` is what one call to a plugin's `compute()` produced for one
replicate, before the framework turns it into a `ReplicateArtifact`: the
reduced observables, the sidecars written for them, and the frame selection
they were measured over.

Both live in `polyzymd.analyses.mda.lifecycle` and are documented on
{doc}`lifecycle`. There is no job object, no function adapter and no backend
policy; a plugin supplies `compute()` and the lifecycle calls it once.
