"""Selection diagnostics."""


def test_a_selection_spanning_a_glycan_segment_suggests_segid_prot(caplog) -> None:
    """A modified protein's glycan shares chain A but is its own segment; a resid-only
    selection that also picks it is flagged, and the fix names the PROT segment."""
    import MDAnalysis as mda
    import numpy as np

    from polyzymd.analyses.shared.diagnostics import warn_if_multi_chain_selection

    universe = mda.Universe.empty(
        2,
        n_residues=2,
        n_segments=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 1],
        trajectory=True,
    )
    universe.add_TopologyAttr("segids", ["PROT", "G001"])
    universe.add_TopologyAttr("resids", [1, 1])
    universe.add_TopologyAttr("resnames", ["ASN", "4YB"])
    universe.add_TopologyAttr("chainIDs", ["A", "A"])
    universe.atoms.positions = np.zeros((2, 3))
    with caplog.at_level("WARNING"):
        assert warn_if_multi_chain_selection(universe.select_atoms("resid 1"), "resid 1")
    assert "segid PROT and" in caplog.text
