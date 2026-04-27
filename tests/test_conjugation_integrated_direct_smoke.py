"""Opt-in direct linkage smoke test without Pablo dependency."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from polyzymd.builders.conjugation.direct_openff import (
    DirectOpenFFLinkageResult,
    build_direct_openff_linkage,
)
from polyzymd.builders.conjugation.linkers import NhsLysModifierLinker
from polyzymd.builders.conjugation.placement import (
    PackmolModifierPlacementSettings,
    place_modifier_with_packmol,
)
from polyzymd.builders.conjugation.polymer_recipe import (
    generate_polymerist_smoke_polymer,
    sbma_nhs_egpma_acb_recipe,
)
from polyzymd.builders.conjugation.polymerist_pdb import generated_fragment_from_polymerist_pdb

SMOKE_ENV_VAR = "POLYZYMD_RUN_CONJUGATION_SMOKE"
POC_PROTEIN_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "polyzymd"
    / "builders"
    / "conjugation"
    / "poc"
    / "data"
    / "NH3_terminal_His_proton_updated.pdb"
)


@pytest.mark.slow
@pytest.mark.conjugation_stack
def test_opt_in_direct_openff_conjugation_smoke(tmp_path: Path):
    """Generate, place, and directly link the ACB polymer to Lys23."""
    _require_direct_smoke_stack_or_skip()
    result = build_direct_smoke_artifact(tmp_path / "direct-conjugation-smoke")

    assert result.linked_pdb_path.exists()
    assert result.summary_json_path.exists()
    assert result.linked_bond.conect_present is True
    assert result.linked_bond.protein_chain_id == "A"
    assert result.linked_bond.modifier_chain_id == "C"
    assert result.finite_coordinates is True


def build_direct_smoke_artifact(artifact_dir: Path) -> DirectOpenFFLinkageResult:
    """Create a direct linked ACB protein-polymer artifact for opt-in smokes."""
    if not POC_PROTEIN_PATH.exists():
        pytest.fail(f"Required lightweight POC protein fixture is missing: {POC_PROTEIN_PATH}")

    recipe = sbma_nhs_egpma_acb_recipe()
    generation = generate_polymerist_smoke_polymer(
        recipe,
        artifact_dir / "polymerist-cache",
        force_regenerate=True,
        max_retries=1,
    )
    if generation.pdb_path is None:
        pytest.fail("Polymerist did not return a generated polymer PDB path")

    modifier = generated_fragment_from_polymerist_pdb(
        generation.pdb_path,
        recipe=recipe,
        sequence=generation.sequence,
        name="direct-acb-smoke",
    )
    linker = NhsLysModifierLinker(target_residue_number=23)
    placement = place_modifier_with_packmol(
        POC_PROTEIN_PATH,
        modifier,
        linker,
        artifact_dir,
        settings=PackmolModifierPlacementSettings(
            tolerance_angstrom=2.0,
            reactive_sphere_radius_angstrom=5.0,
            nloop=500,
            movebadrandom=True,
            target_bond_length_angstrom=1.33,
        ),
    )
    plan = linker.resolve_plan(POC_PROTEIN_PATH, modifier)
    return build_direct_openff_linkage(
        protein_pdb_path=POC_PROTEIN_PATH,
        modifier=placement.placed_modifier,
        resolved_plan=plan,
        output_dir=artifact_dir,
    )


def _require_direct_smoke_stack_or_skip() -> None:
    """Skip unless the opt-in direct construction stack is available."""
    if os.environ.get(SMOKE_ENV_VAR) != "1":
        pytest.skip(f"Set {SMOKE_ENV_VAR}=1 to run the direct conjugation smoke")
    if shutil.which("packmol") is None:
        pytest.skip("Packmol binary is not available on PATH")

    from polyzymd.builders.conjugation.polymerist_compat import ensure_polymerist_py312_compat

    ensure_polymerist_py312_compat()
    pytest.importorskip("polymerist", exc_type=ImportError)
    pytest.importorskip("openff.toolkit")
    pytest.importorskip("rdkit")
