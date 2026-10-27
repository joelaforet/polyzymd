"""dynamics_seed: one seed of the dynamics per replicate and phase."""

from __future__ import annotations


def test_a_seed_depends_on_the_replicate_and_the_phase_only() -> None:
    """The seed is reproducible, differs between replicates and phases, and lies in range."""
    from polyzymd.simulation.seeds import MAX_SEED, dynamics_seed

    assert dynamics_seed(1, "production:0") == dynamics_seed(1, "production:0")
    seeds = {
        dynamics_seed(r, p)
        for r in (1, 2, 3)
        for p in ("velocities", "production:0", "production:1")
    }
    assert len(seeds) == 9
    assert all(1 <= seed <= MAX_SEED for seed in seeds)
