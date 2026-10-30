"""Mutate, protonate and covalently modify proteins, recording every change.

Usable on its own, without a PolyzyMD config::

    from polyzymd import modifying_proteins as mp

    protein = mp.load_protein("ubiquitin.pdb")
    mp.deprotonate(protein, "A:K63", "NZ")
    acyl = mp.fragment_from_smiles("*C(=O)CCCCCCC", "OC8")
    mp.attach(protein, acyl, type="generic", site="A:K63", atom_name="NZ")

The names and arguments are PolyzyMD's and stay stable while the mBuild
biopolymers API they call changes.
"""

from polyzymd.modifying_proteins.api import (
    ATTACHMENT_TYPES,
    AtomKey,
    ModificationPlan,
    ModificationRecord,
    ModifiedProtein,
    apply,
    attach,
    deprotonate,
    fragment_from_ccd,
    fragment_from_pdb,
    fragment_from_smiles,
    fragment_source,
    load_protein,
    mutate,
    plan_attach,
    plan_deprotonate,
    plan_mutate,
    plan_protonate,
    protonate,
    relax_all,
)
from polyzymd.modifying_proteins.parameters import (
    ParameterizationRecord,
    ParameterSet,
    SiteParameters,
    parameterize,
)
from polyzymd.modifying_proteins.sites import Site, expand_sites, parse_site
from polyzymd.modifying_proteins.spec import apply_spec, expand, to_spec, to_yaml
from polyzymd.modifying_proteins.topology import ModificationSite, modification_sites, segments

__all__ = [
    "ATTACHMENT_TYPES",
    "AtomKey",
    "ModificationPlan",
    "ModificationRecord",
    "ModificationSite",
    "ModifiedProtein",
    "ParameterSet",
    "ParameterizationRecord",
    "Site",
    "SiteParameters",
    "apply",
    "apply_spec",
    "attach",
    "deprotonate",
    "expand",
    "expand_sites",
    "fragment_from_ccd",
    "fragment_from_pdb",
    "fragment_from_smiles",
    "fragment_source",
    "load_protein",
    "modification_sites",
    "mutate",
    "parameterize",
    "parse_site",
    "plan_attach",
    "plan_deprotonate",
    "plan_mutate",
    "plan_protonate",
    "protonate",
    "relax_all",
    "segments",
    "to_spec",
    "to_yaml",
]
