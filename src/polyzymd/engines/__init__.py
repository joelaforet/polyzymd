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


def create_engine(
    config: object,
    override: str | None = None,
    defer_binary: bool = False,
) -> SimulationEngine:
    """Create an engine instance from configuration.

    Parameters
    ----------
    config : object
        Simulation configuration object.
    override : str | None, optional
        CLI-provided engine override.
    defer_binary : bool, optional
        When True, allow engine-specific deferred binary resolution for
        scheduler-only paths.

    Returns
    -------
    SimulationEngine
        Instantiated engine backend.

    Raises
    ------
    ValueError
        If neither an explicit override nor a valid ``config.engine`` value is
        provided.
    """
    if override is not None:
        if not isinstance(override, str) or not override.strip():
            raise ValueError("Engine override must be a non-empty string")
        engine_name = override.strip()
    else:
        configured_engine = getattr(config, "engine", None)
        if not isinstance(configured_engine, str) or not configured_engine.strip():
            raise ValueError(
                "Simulation config must define a non-empty string engine "
                "('openmm' or 'gromacs'), or pass an explicit engine override"
            )
        engine_name = configured_engine.strip()
    cls = get_engine_class(engine_name)
    if engine_name.lower() == "gromacs":
        return cls.from_config(config, defer_binary=defer_binary)
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
