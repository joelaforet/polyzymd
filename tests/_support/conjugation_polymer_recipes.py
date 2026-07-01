"""Test support recipes for conjugation polymer fixtures."""

from __future__ import annotations

from polyzymd.builders.conjugation.polymer.recipe import PolymerMonomerRecipe, PolymerRecipe

SBMA_SMILES = (
    "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+]"
    "(C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])"
    "C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
)
EGPMA_SMILES = (
    "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
)
NHS_SMILES = "CC(=C)C(=O)ON1C(=O)CCC1=O"


def sbma_egpma_nhs_recipe(
    *,
    length: int = 10,
    seed: int | None = 42,
    reactive_monomer_index: int | None = None,
    fixed_sequence: str | None = None,
) -> PolymerRecipe:
    """Build an SBMA/EGPMA/NHS recipe fixture.

    Parameters
    ----------
    length : int, optional
        Degree of polymerization, by default 10.
    seed : int or None, optional
        Random seed for deterministic sequence generation, by default 42.
    reactive_monomer_index : int or None, optional
        Zero-based NHS residue index. ``None`` centers NHS, by default ``None``.
    fixed_sequence : str or None, optional
        Exact monomer-label sequence overriding stochastic generation, by default
        ``None``.

    Returns
    -------
    PolymerRecipe
        Validated stochastic polymer recipe.
    """
    return PolymerRecipe(
        name="SBMA-EGPMA-NHS",
        monomers=(
            PolymerMonomerRecipe(
                label="A",
                name="SBMA",
                residue_name="SBM",
                smiles=SBMA_SMILES,
                probability=0.945,
            ),
            PolymerMonomerRecipe(
                label="B",
                name="EGPMA",
                residue_name="EGP",
                smiles=EGPMA_SMILES,
                probability=0.045,
            ),
            PolymerMonomerRecipe(
                label="C",
                name="NHS",
                residue_name="NHS",
                smiles=NHS_SMILES,
                probability=0.01,
            ),
        ),
        length=length,
        seed=seed,
        reactive_monomer_label="C",
        reactive_monomer_index=reactive_monomer_index,
        fixed_sequence=fixed_sequence,
    )


def sbma_nhs_egpma_acb_recipe() -> PolymerRecipe:
    """Build the deterministic SBMA:NHS:EGPMA recipe fixture.

    Returns
    -------
    PolymerRecipe
        Three-monomer recipe whose fixed sequence ``ACB`` maps to
        SBMA:NHS:EGPMA with the NHS monomer centered for Lys linkage.
    """
    return sbma_egpma_nhs_recipe(
        length=3,
        seed=None,
        reactive_monomer_index=1,
        fixed_sequence="ACB",
    )
