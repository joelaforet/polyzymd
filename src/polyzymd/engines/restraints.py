"""Engine-neutral restraint specifications for simulation engines."""

from pydantic import BaseModel, Field


class RestraintSpec(BaseModel):
    """Base restraint specification.

    Parameters
    ----------
    group : str
        Logical atom group name the restraint applies to.
    force_constant : float
        Restraint force constant.
    force_constant_unit : str, optional
        Force-constant unit string, by default ``"kJ/mol/nm^2"``.
    """

    group: str
    force_constant: float = Field(gt=0.0)
    force_constant_unit: str = "kJ/mol/nm^2"


class PositionRestraintSpec(RestraintSpec):
    """Position restraint specification.

    Parameters
    ----------
    reference_positions : str, optional
        Reference coordinate source. Use ``"initial"`` or a path string.
    """

    reference_positions: str = "initial"


class FlatBottomRestraintSpec(RestraintSpec):
    """Flat-bottom restraint specification placeholder.

    Parameters
    ----------
    radius : float
        Flat-bottom radius value.
    radius_unit : str, optional
        Radius unit string, by default ``"nanometer"``.
    """

    radius: float = Field(gt=0.0)
    radius_unit: str = "nanometer"
