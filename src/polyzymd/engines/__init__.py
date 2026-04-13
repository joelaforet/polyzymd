"""Multi-engine simulation dispatch for PolyzyMD."""

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine, TrajectoryLayout


def get_engine_class(name: str) -> type[SimulationEngine]:
    """Look up an engine class by name.

    Parameters
    ----------
    name : str
        Engine name.

    Returns
    -------
    type[SimulationEngine]
        Engine implementation class.

    Raises
    ------
    ValueError
        If the engine name is not recognized.
    """
    match name.lower():
        case "openmm":
            from polyzymd.engines.openmm import OpenMMEngine

            return OpenMMEngine
        case "gromacs":
            from polyzymd.engines.gromacs import GromacsEngine

            return GromacsEngine
        case _:
            raise ValueError(f"Unknown engine: {name!r}. Available: openmm, gromacs")


def create_engine(config: object, override: str | None = None) -> SimulationEngine:
    """Create an engine instance from configuration.

    Parameters
    ----------
    config : object
        Simulation configuration object.
    override : str | None, optional
        CLI-provided engine override.

    Returns
    -------
    SimulationEngine
        Instantiated engine backend.
    """
    engine_name = override or getattr(config, "engine", "openmm") or "openmm"
    cls = get_engine_class(engine_name)
    return cls.from_config(config)


def list_engines() -> dict[str, str]:
    """Return available engine names and descriptions.

    Returns
    -------
    dict[str, str]
        Mapping of engine names to human-readable descriptions.
    """
    return {
        "openmm": "OpenMM GPU-accelerated molecular dynamics",
        "gromacs": "GROMACS molecular dynamics",
    }


__all__ = [
    "SimulationEngine",
    "TrajectoryLayout",
    "EngineSubmitRequest",
    "get_engine_class",
    "create_engine",
    "list_engines",
]
