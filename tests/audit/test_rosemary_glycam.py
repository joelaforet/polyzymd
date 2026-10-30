"""The Rosemary/GLYCAM06j hybrid: GLYCAM06j glycans on a protein parameterized whole by Rosemary.

An expert option for benchmarking (design decision of 2026-10-30): each glycan
site takes GLYCAM06j's charges, terms and 1-4 pairs, and the charge that swap
removes or adds is spread over the linkage residue's backbone.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("mbuild.biopolymers")

from polyzymd.audit.findings import run_checks  # noqa: E402
from polyzymd.audit.realized import RealizedSystem  # noqa: E402
from polyzymd.modifying_proteins import ParameterSet  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.rosemary import newest_rosemary  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
CHITOBIOSE = ROOT / "tests" / "data" / "modifying_proteins" / "glycam_G42666HT.pdb"
GLYCAM = ParameterSet(["GLYCAM_06j-1.xml"], {"method": "library"})
ROSEMARY = newest_rosemary()
NAGL = {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"}
#: The O-glycan on Ser13 needs about 0.0095 e per backbone atom.
BOUND = 0.01
CHECKS = [
    "charges.source",
    "charges.site_integer",
    "charges.site_shift",
    "parameters.bonded",
    "parameters.nonbonded",
    "parameters.one_four",
    "parameters.modified_backbone",
    "parameters.glycan_force_field",
    "masses",
    "modifications.sites",
]

pytestmark = pytest.mark.filterwarnings("ignore")


def _glycoprotein():
    """Trp-cage with an N-glycan on Asn5 (from Q5N) and an O-glycan on Ser13."""
    from polyzymd import modifying_proteins as mp

    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:Q5", "ASN")
    for site, atom_name, reaction in (
        ("A:N5", "ND2", "N-glycosylation"),
        ("A:S13", "OG", "O-glycosylation"),
    ):
        mp.attach(
            protein,
            {"pdb_path": str(CHITOBIOSE)},
            type="glycan",
            site=site,
            atom_name=atom_name,
            reaction=reaction,
            fragment_atom_name="C1",
            fragment_resnum=2,
        )
    return protein


@pytest.fixture(scope="module")
def hybrid():
    from polyzymd import modifying_proteins as mp

    protein = _glycoprotein()
    result = mp.parameterize_whole(
        protein, (ROSEMARY,), NAGL, glycan=GLYCAM, max_charge_shift=BOUND
    )
    return protein, result


def _realize(protein, result, folder: Path, edit=None):
    import openmm

    from polyzymd.modifying_proteins.components import protein_component

    system = result.to_openmm()
    if edit is not None:
        edit(system)
    (folder / "system.xml").write_text(openmm.XmlSerializer.serialize(system))
    protein.save_pdb(folder / "protein.pdb")
    realized = RealizedSystem.from_openmm(folder / "system.xml", folder / "protein.pdb")
    component = protein_component(protein, result.record, {"glycan": GLYCAM}, name="glyco")
    return realized, component


def test_glycan_sites_take_glycam_and_the_backbone_keeps_the_charge(hybrid) -> None:
    protein, result = hybrid
    record = result.record
    assert [site.name for site in record.sites] == ["G001", "G002"]
    assert [site.backbone_residues for site in record.sites] == [["ASN"], ["SER"]]
    assert len(record.glycam) == 2 and all(g.n_pairs > 0 for g in record.glycam)
    assert record.total_charge == pytest.approx(1.0, abs=1e-6)
    for site in record.sites:
        # The shift moves NAGL's site sum minus GLYCAM's onto six backbone atoms.
        assert site.shift_total == pytest.approx(site.nagl_sum - site.target)
        assert site.shift_per_atom == pytest.approx(site.shift_total / 6)
        assert abs(site.shift_per_atom) <= BOUND


def test_a_backbone_shift_above_the_bound_is_refused() -> None:
    from polyzymd import modifying_proteins as mp

    with pytest.raises(ValueError, match="Raise max_charge_shift"):
        mp.parameterize_whole(_glycoprotein(), (ROSEMARY,), NAGL, glycan=GLYCAM)


def test_the_audit_passes_the_hybrid_and_warns_that_it_is_one(hybrid, tmp_path) -> None:
    realized, component = _realize(*hybrid, tmp_path)
    findings = {f.check: f for f in run_checks(realized, [component], CHECKS)}
    not_passing = {
        check: (f.status, f.detail)
        for check, f in findings.items()
        if f.status not in ("PASS", "INFO") and check != "parameters.glycan_force_field"
    }
    assert not_passing == {}
    warning = findings["parameters.glycan_force_field"]
    assert warning.status == "WARN" and "Rosemary/GLYCAM06j hybrid" in warning.detail


def test_the_audit_fails_an_unshifted_linkage_backbone(hybrid, tmp_path) -> None:
    import openmm

    protein, result = hybrid
    serine = result.record.sites[1]
    keys = protein.atom_keys()
    linkage = [
        i
        for i, key in enumerate(keys)
        if (key.chain_id, key.resnum) == ("A", 13)
        and key.atom_name in ("N", "CA", "C")
        and i not in serine.atom_indices
    ]

    def undo_shift(system):
        force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
        for i in linkage:
            q, sigma, epsilon = force.getParticleParameters(i)
            force.setParticleParameters(i, q._value - serine.shift_per_atom, sigma, epsilon)

    realized, component = _realize(protein, result, tmp_path, edit=undo_shift)
    (finding,) = run_checks(realized, [component], ["parameters.modified_backbone"])
    assert finding.status == "FAIL"
    assert np.isclose(abs(serine.shift_per_atom), 0.0095, atol=0.002)
