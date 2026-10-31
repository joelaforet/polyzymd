"""Sizing of periodic boxes, and the rotation that precedes a rectangular one."""

import numpy as np
import pytest

from polyzymd.utils.boxvectors import (
    describe_box_plan,
    plan_box,
    plan_rectangular_box,
    principal_axes_rotation,
)


def _rod(length_nm: float = 8.0, n: int = 400, seed: int = 3) -> np.ndarray:
    """A flattened rod (axes 0.3, 1.0 and ``length_nm`` wide), turned to an arbitrary angle."""
    rng = np.random.default_rng(seed)
    body = rng.uniform(-0.5, 0.5, size=(n, 3)) * np.array([0.3, 1.0, length_nm])
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1.0
    return body @ q.T + np.array([1.0, -2.0, 3.0])


def test_the_longest_principal_axis_ends_up_on_z_and_distances_are_kept() -> None:
    points = _rod()
    rotation = principal_axes_rotation(points)
    assert np.isclose(np.linalg.det(rotation), 1.0)
    center = points.mean(axis=0)
    turned = (points - center) @ rotation.T + center
    extent = turned.max(axis=0) - turned.min(axis=0)
    assert extent[2] == pytest.approx(extent.max())
    assert extent[0] == pytest.approx(extent.min())
    before = np.linalg.norm(points[:, None] - points[None], axis=-1)
    after = np.linalg.norm(turned[:, None] - turned[None], axis=-1)
    assert np.allclose(before, after)


def test_a_rectangular_box_is_the_bounding_box_plus_twice_each_padding() -> None:
    points = np.array([[0.0, 0.0, 0.0], [2.0, 1.0, 8.0]])
    plan = plan_rectangular_box(points, (1.0, 1.5, 2.0), margin_nm=0.2)
    assert np.allclose(np.diagonal(plan["box_vectors"]), [4.0, 4.0, 12.0])
    assert np.allclose(plan["box_vectors"], np.diag(np.diagonal(plan["box_vectors"])))
    assert "rectangular 4.00 x 4.00 x 12.00 nm" in describe_box_plan(plan)


def test_a_long_solute_needs_far_less_volume_in_an_aligned_rectangular_box() -> None:
    """Why the option exists: a dodecahedron is sized by the solute's longest dimension."""
    import openff.packmol as packmol

    points = _rod(length_nm=15.0)
    center = points.mean(axis=0)
    aligned = (points - center) @ principal_axes_rotation(points).T + center
    rectangular = plan_rectangular_box(aligned, (1.0, 1.0, 1.0), margin_nm=0.2)
    dodecahedron = plan_box(points, packmol.RHOMBIC_DODECAHEDRON, padding_nm=1.0, margin_nm=0.2)
    volume = lambda plan: abs(np.linalg.det(plan["box_vectors"]))  # noqa: E731
    assert volume(rectangular) < 0.35 * volume(dodecahedron)


def _rows_have_canonical_signs(rotation: np.ndarray) -> bool:
    return all(row[np.argmax(np.abs(row))] > 0 for row in rotation[1:])


@pytest.mark.parametrize("seed", range(8))
def test_the_rotation_has_canonical_signs_and_ignores_atom_order_and_noise(seed: int) -> None:
    points = _rod(seed=seed)
    rotation = principal_axes_rotation(points)
    assert np.isclose(np.linalg.det(rotation), 1.0)
    # y and z keep the sign rule; x is the axis flipped, if needed, to stay proper.
    assert _rows_have_canonical_signs(rotation)
    rng = np.random.default_rng(seed + 100)
    shuffled = points[rng.permutation(len(points))]
    noisy = points + rng.normal(scale=1e-6, size=points.shape)
    assert np.allclose(principal_axes_rotation(shuffled), rotation, atol=1e-6)
    assert np.allclose(principal_axes_rotation(noisy), rotation, atol=1e-4)


def _disk(tilt: float) -> np.ndarray:
    """A thin disk in the input y-z plane, its normal tilted by ``tilt`` rad: its two
    smallest moments are within 1% of each other."""
    radius, angle = np.meshgrid(np.linspace(0.1, 3.0, 30), np.linspace(0.0, 2.0 * np.pi, 72, False))
    radius, angle = radius.ravel(), angle.ravel() + 0.003  # not symmetric about y or z
    disk = np.column_stack([np.zeros_like(radius), radius * np.cos(angle), radius * np.sin(angle)])
    disk = np.vstack([disk + [0.1, 0.0, 0.0], disk - [0.1, 0.0, 0.0]])
    turn = np.array(
        [[np.cos(tilt), -np.sin(tilt), 0.0], [np.sin(tilt), np.cos(tilt), 0.0], [0, 0, 1]]
    )
    return disk @ turn.T


def test_a_degenerate_plane_keeps_the_input_orientation(caplog) -> None:
    import logging

    with caplog.at_level(logging.INFO):
        rotation = principal_axes_rotation(_disk(0.0))
    assert "within 1%" in caplog.text
    # x is the disk normal; within the disk, y and z stay where they were.
    assert np.allclose(rotation, np.eye(3), atol=1e-6)
    # A slightly tilted disk turns only by the tilt, not about its normal.
    tilted = principal_axes_rotation(_disk(0.01))
    assert np.allclose(tilted, rotation, atol=0.02)


def test_grafted_chains_do_not_turn_the_box() -> None:
    """Replicates draw their own grafted chains; the rotation must not depend on them."""
    from polyzymd.utils.boxvectors import fixed_solute_alignment

    protein = _rod(seed=5)
    rng = np.random.default_rng(9)
    rotations = []
    for _ in range(2):
        chain = rng.uniform(-3.0, 3.0, size=(60, 3)) + np.array([4.0, 0.0, 0.0])
        positions = np.vstack([protein, chain])
        exclude = np.r_[np.zeros(len(protein), bool), np.ones(len(chain), bool)]
        rotation, center = fixed_solute_alignment(
            positions, np.full(len(positions), 12.0), np.full(len(positions), 6), exclude
        )
        rotations.append(rotation)
        assert np.allclose(center, protein.mean(axis=0))
    assert np.allclose(rotations[0], rotations[1])
