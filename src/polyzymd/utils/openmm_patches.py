"""
Monkey-patches for OpenFF Interchange performance bugs.

This module contains targeted fixes for known performance issues in
OpenFF Interchange that have not yet been fixed upstream. Each patch is
version-guarded and idempotent.

Patches
-------
patch_interchange_to_openmm_14_pairs
    Fixes O(N^2) list membership check in ``_create_multiple_nonbonded_forces``.
    Affects: Interchange 0.5.x
    Upstream: https://github.com/openmm/openff-interchange (pending issue)

Notes
-----
These patches are applied at runtime before calling ``to_openmm()``. They
modify module-level function references, so they affect all subsequent calls
within the process. The original functions are preserved and can be restored.
"""

import logging
from typing import Any

logger = logging.getLogger(__name__)

# Sentinel to track whether the patch has been applied
_PATCH_APPLIED = False
_ORIGINAL_FN: Any = None


def _make_patched_create_multiple_nonbonded_forces():
    """Build the patched version of _create_multiple_nonbonded_forces.

    The ONLY semantic change from Interchange 0.5.0:
        Line 682: ``list()`` -> ``set()``
        Line 695: ``.append()`` -> ``.add()``

    This changes the 1-4 pair membership check from O(N) to O(1),
    reducing the overall complexity of the exception loop from
    O(exceptions * 1_4_pairs) to O(exceptions).

    Returns
    -------
    callable
        The patched function with the same signature as the original.
    """
    import openmm
    import openmm.unit  # noqa: F401
    from openff.interchange.components.toolkit import _get_14_pairs
    from openff.interchange.exceptions import (
        UnsupportedCutoffMethodError,
        UnsupportedExportError,
    )
    from openff.interchange.interop.openmm._nonbonded import (
        _create_electrostatics_force,
        _create_vdw_force,
        _get_14_scaling_factors,
        _get_scaled_potential_function,
        _set_particle_parameters,
    )

    def _create_multiple_nonbonded_forces_patched(
        data,
        interchange,
        system,
        ewald_tolerance,
        molecule_virtual_site_map,
        openff_openmm_particle_map,
    ):
        if molecule_virtual_site_map in (None, {}):
            has_virtual_sites = False
        elif all(len(v) == 0 for v in molecule_virtual_site_map.values()):
            has_virtual_sites = False
        else:
            has_virtual_sites = True

        vdw_force = _create_vdw_force(
            data,
            interchange,
            molecule_virtual_site_map,
            has_virtual_sites,
        )

        electrostatics_force = _create_electrostatics_force(
            data,
            interchange,
            ewald_tolerance,
            molecule_virtual_site_map,
            has_virtual_sites,
            openff_openmm_particle_map,
        )

        _set_particle_parameters(
            data,
            vdw_force,
            electrostatics_force,
            interchange,
            has_virtual_sites,
            molecule_virtual_site_map,
            openff_openmm_particle_map,
        )

        coul_const = 138.935456  # kJ/nm

        if vdw_force is not None:
            vdw = data.vdw_collection

            if vdw.is_plugin:
                vdw_14_force = openmm.CustomBondForce(
                    _get_scaled_potential_function(data.vdw_expression),
                )
                vdw_14_force.setName("vdW 1-4 force")
                for index in [1, 2]:
                    for parameter in vdw.potential_parameters():
                        vdw_14_force.addPerBondParameter(f"{parameter}{index!s}")
                vdw_14_force.addGlobalParameter("scale14", vdw.scale_14)
                for global_parameter in vdw.global_parameters():
                    vdw_14_force.addGlobalParameter(
                        global_parameter,
                        getattr(vdw, global_parameter).m,
                    )
                for term, value in data.vdw_collection.pre_computed_terms().items():
                    vdw_14_force.addGlobalParameter(term, value)
            else:
                vdw_expression: str = data.vdw_expression
                vdw_14_force = openmm.CustomBondForce(vdw_expression)
                vdw_14_force.setName("vdW 1-4 force")
                for parameter in vdw.potential_parameters():
                    vdw_14_force.addPerBondParameter(f"{parameter}")

            vdw_14_force.setUsesPeriodicBoundaryConditions(interchange.box is not None)
        else:
            vdw_14_force = None

        coul_14_force = openmm.CustomBondForce(f"{coul_const}*qq/r")
        coul_14_force.setName("Electrostatics 1-4 force")
        coul_14_force.addPerBondParameter("qq")
        coul_14_force.setUsesPeriodicBoundaryConditions(interchange.box is not None)

        coul_14, vdw_14 = _get_14_scaling_factors(data)

        # ==================================================================
        # PATCHED: set() instead of list() for O(1) membership checks.
        # Original Interchange 0.5.0 uses list() + .append(), causing O(N)
        # linear scans in the exception loop below. For large systems
        # (>100K atoms), this single change reduces to_openmm() from
        # ~53 min to ~4 min.
        # ==================================================================
        openmm_pairs = set()

        for atom1, atom2 in _get_14_pairs(interchange.topology):
            openff_indices = (
                interchange.topology.atom_index(atom1),
                interchange.topology.atom_index(atom2),
            )
            openmm_indices = (
                openff_openmm_particle_map[openff_indices[0]],
                openff_openmm_particle_map[openff_indices[1]],
            )
            openmm_pairs.add(openmm_indices)

        if electrostatics_force is not None:
            for i in range(electrostatics_force.getNumExceptions()):
                (p1, p2, _, _, _) = electrostatics_force.getExceptionParameters(i)

                if (p1, p2) in openmm_pairs or (p2, p1) in openmm_pairs:
                    if vdw_force is not None:
                        if data.vdw_collection.is_plugin:
                            parameters1 = vdw_force.getParticleParameters(p1)
                            parameters2 = vdw_force.getParticleParameters(p2)
                            vdw_14_force.addBond(p1, p2, [*parameters1, *parameters2])
                        else:
                            sig1, eps1 = vdw_force.getParticleParameters(p1)
                            sig2, eps2 = vdw_force.getParticleParameters(p2)

                            if data.mixing_rule == "lorentz-berthelot":
                                sig_14 = (sig1 + sig2) * 0.5
                                eps_14 = (eps1 * eps2) ** 0.5 * vdw_14
                            elif data.mixing_rule == "geometric":
                                sig_14 = (sig1 * sig2) ** 0.5
                                eps_14 = (eps1 * eps2) ** 0.5 * vdw_14
                            else:
                                raise UnsupportedExportError(
                                    f"Unsupported mixing rule: {data.mixing_rule}",
                                )

                            vdw_14_force.addBond(p1, p2, [sig_14, eps_14])

                    q1 = electrostatics_force.getParticleParameters(p1)[0]
                    q2 = electrostatics_force.getParticleParameters(p2)[0]
                    qq = q1 * q2 * coul_14
                    coul_14_force.addBond(p1, p2, [qq])

                if vdw_force is not None:
                    vdw_force.addExclusion(p1, p2)

                if electrostatics_force is not None:
                    electrostatics_force.setExceptionParameters(i, p1, p2, 0.0, 0.0, 0.0)

        for force in [vdw_force, electrostatics_force, vdw_14_force, coul_14_force]:
            if force is not None:
                system.addForce(force)

        if vdw_force is not None and electrostatics_force is not None:
            vdw_uses = vdw_force.usesPeriodicBoundaryConditions()
            elec_uses = electrostatics_force.usesPeriodicBoundaryConditions()
            if vdw_uses != elec_uses:
                raise UnsupportedCutoffMethodError(
                    "When using `openmm.CustomNonbondedForce`, vdW and "
                    "electrostatics cutoff methods must agree on whether or "
                    "not periodic boundary conditions should be used. Found "
                    f"vdw method {vdw_force.getNonbondedMethod()}, and "
                    f"electrostatics method "
                    f"{electrostatics_force.getNonbondedMethod()}, ",
                )

    return _create_multiple_nonbonded_forces_patched


def patch_interchange_to_openmm_14_pairs() -> bool:
    """Apply the list-to-set patch for 1-4 pair lookup in to_openmm().

    This patches ``_create_multiple_nonbonded_forces`` in the OpenFF
    Interchange ``openmm._nonbonded`` module to use a ``set`` instead of
    a ``list`` for storing 1-4 pair indices. This changes the membership
    check from O(N) to O(1), fixing an O(N^2) bottleneck in
    ``to_openmm(combine_nonbonded_forces=False)``.

    The patch is:
    - **Version-guarded**: Only applied to Interchange 0.5.x
    - **Idempotent**: Safe to call multiple times; only patches once
    - **Reversible**: Call ``unpatch_interchange_to_openmm_14_pairs()``

    Returns
    -------
    bool
        True if the patch was applied (or was already applied), False if
        the Interchange version is not 0.5.x.
    """
    global _PATCH_APPLIED, _ORIGINAL_FN

    if _PATCH_APPLIED:
        return True

    try:
        import openff.interchange

        version = openff.interchange.__version__
    except (ImportError, AttributeError):
        logger.debug("OpenFF Interchange not installed; skipping to_openmm patch")
        return False

    if not version.startswith("0.5"):
        logger.info(
            "Interchange version %s is not 0.5.x; skipping to_openmm "
            "list-to-set patch (may be fixed upstream)",
            version,
        )
        return False

    try:
        import openff.interchange.interop.openmm._nonbonded as nb_module

        _ORIGINAL_FN = nb_module._create_multiple_nonbonded_forces
        patched_fn = _make_patched_create_multiple_nonbonded_forces()
        nb_module._create_multiple_nonbonded_forces = patched_fn
        _PATCH_APPLIED = True

        logger.info(
            "Applied to_openmm list-to-set patch for Interchange %s "
            "(1-4 pair lookup: O(N) -> O(1))",
            version,
        )
        return True
    except Exception:
        logger.warning(
            "Failed to apply to_openmm list-to-set patch",
            exc_info=True,
        )
        return False


def unpatch_interchange_to_openmm_14_pairs() -> bool:
    """Restore the original _create_multiple_nonbonded_forces function.

    Returns
    -------
    bool
        True if the original was restored, False if no patch was active.
    """
    global _PATCH_APPLIED, _ORIGINAL_FN

    if not _PATCH_APPLIED or _ORIGINAL_FN is None:
        return False

    try:
        import openff.interchange.interop.openmm._nonbonded as nb_module

        nb_module._create_multiple_nonbonded_forces = _ORIGINAL_FN
        _PATCH_APPLIED = False
        _ORIGINAL_FN = None
        logger.info("Restored original _create_multiple_nonbonded_forces")
        return True
    except Exception:
        logger.warning(
            "Failed to restore original _create_multiple_nonbonded_forces",
            exc_info=True,
        )
        return False
