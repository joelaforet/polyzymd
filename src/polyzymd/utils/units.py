"""
Unit conversion utilities between OpenFF and OpenMM unit systems.

This module provides helpers for converting quantities between:
- OpenFF: Uses openff.units.Quantity (based on Pint)
- OpenMM: Uses openmm.unit.Quantity (OpenMM's built-in unit system)

These conversions are necessary when interfacing between OpenFF Interchange
objects (which use OpenFF units) and OpenMM simulation objects.

Example:
    >>> from openff.units import Quantity as OFFQuantity
    >>> from polyzymd.utils.units import openff_to_openmm, openmm_to_openff
    >>>
    >>> # Convert OpenFF quantity to OpenMM
    >>> off_distance = OFFQuantity(1.5, "nanometer")
    >>> omm_distance = openff_to_openmm(off_distance)
    >>>
    >>> # Convert OpenMM quantity to OpenFF
    >>> import openmm.unit as omm_unit
    >>> omm_energy = 100.0 * omm_unit.kilojoule_per_mole
    >>> off_energy = openmm_to_openff(omm_energy)
"""

from typing import Union

# Type aliases for clarity
# These will be actual types at runtime when the packages are installed
OpenFFQuantity = "openff.units.Quantity"
OpenMMQuantity = "openmm.unit.Quantity"


def openff_to_openmm(quantity: "OpenFFQuantity") -> "OpenMMQuantity":
    """Convert an OpenFF Quantity to an OpenMM Quantity.

    This function handles the conversion between Pint-based OpenFF units
    and OpenMM's built-in unit system. It extracts the magnitude and unit
    string from the OpenFF quantity and reconstructs it as an OpenMM quantity.

    Args:
        quantity: An OpenFF Quantity (openff.units.Quantity) with units.

    Returns:
        The equivalent OpenMM Quantity (openmm.unit.Quantity).

    Raises:
        ValueError: If the unit cannot be converted.

    Example:
        >>> from openff.units import Quantity
        >>> off_qty = Quantity(1.5, "nanometer")
        >>> omm_qty = openff_to_openmm(off_qty)
        >>> print(omm_qty)  # 1.5 nm
    """
    import openmm.unit as omm_unit

    # Get the magnitude (dimensionless value)
    magnitude = quantity.magnitude

    # Get the unit string and convert to OpenMM unit
    # OpenFF uses Pint unit strings like "nanometer", "kilojoule/mole", etc.
    unit_str = str(quantity.units)

    # Map common OpenFF unit strings to OpenMM units
    omm_unit_obj = _get_openmm_unit(unit_str)

    return magnitude * omm_unit_obj


def openmm_to_openff(quantity: "OpenMMQuantity") -> "OpenFFQuantity":
    """Convert an OpenMM Quantity to an OpenFF Quantity.

    This function handles the conversion from OpenMM's built-in unit system
    to Pint-based OpenFF units. It extracts the magnitude and unit from the
    OpenMM quantity and reconstructs it as an OpenFF quantity.

    Args:
        quantity: An OpenMM Quantity (openmm.unit.Quantity) with units.

    Returns:
        The equivalent OpenFF Quantity (openff.units.Quantity).

    Raises:
        ValueError: If the unit cannot be converted.

    Example:
        >>> import openmm.unit as omm_unit
        >>> omm_qty = 1.5 * omm_unit.nanometer
        >>> off_qty = openmm_to_openff(omm_qty)
        >>> print(off_qty)  # 1.5 nanometer
    """
    from openff.units import Quantity as OFFQuantity

    # Get the value in the quantity's own units
    magnitude = quantity.value_in_unit(quantity.unit)

    # Get the unit string from OpenMM
    unit_str = str(quantity.unit)

    # Map OpenMM unit string to OpenFF unit string
    off_unit_str = _get_openff_unit_string(unit_str)

    return OFFQuantity(magnitude, off_unit_str)


def _get_openmm_unit(unit_str: str) -> "openmm.unit.Unit":
    """Map an OpenFF/Pint unit string to an OpenMM unit object.

    Args:
        unit_str: Unit string from OpenFF (e.g., "nanometer", "kilojoule/mole").

    Returns:
        The corresponding OpenMM unit object.

    Raises:
        ValueError: If the unit string is not recognized.
    """
    import openmm.unit as omm_unit

    # Normalize the unit string (lowercase, strip whitespace)
    unit_str = unit_str.lower().strip()

    # Common unit mappings
    unit_map = {
        # Length
        "nanometer": omm_unit.nanometer,
        "nanometers": omm_unit.nanometer,
        "nm": omm_unit.nanometer,
        "angstrom": omm_unit.angstrom,
        "angstroms": omm_unit.angstrom,
        "å": omm_unit.angstrom,
        "picometer": omm_unit.picometer,
        "picometers": omm_unit.picometer,
        "meter": omm_unit.meter,
        "meters": omm_unit.meter,
        # Time
        "picosecond": omm_unit.picosecond,
        "picoseconds": omm_unit.picosecond,
        "ps": omm_unit.picosecond,
        "femtosecond": omm_unit.femtosecond,
        "femtoseconds": omm_unit.femtosecond,
        "fs": omm_unit.femtosecond,
        "nanosecond": omm_unit.nanosecond,
        "nanoseconds": omm_unit.nanosecond,
        "ns": omm_unit.nanosecond,
        # Energy
        "kilojoule / mole": omm_unit.kilojoule_per_mole,
        "kilojoule/mole": omm_unit.kilojoule_per_mole,
        "kj/mol": omm_unit.kilojoule_per_mole,
        "kilojoule_per_mole": omm_unit.kilojoule_per_mole,
        "kilocalorie / mole": omm_unit.kilocalorie_per_mole,
        "kilocalorie/mole": omm_unit.kilocalorie_per_mole,
        "kcal/mol": omm_unit.kilocalorie_per_mole,
        "kilocalorie_per_mole": omm_unit.kilocalorie_per_mole,
        # Force
        "kilojoule / mole / nanometer": omm_unit.kilojoule_per_mole / omm_unit.nanometer,
        "kilojoule / (mole * nanometer)": omm_unit.kilojoule_per_mole / omm_unit.nanometer,
        "kj/mol/nm": omm_unit.kilojoule_per_mole / omm_unit.nanometer,
        # Spring constant (for restraints)
        "kilojoule / mole / nanometer ** 2": omm_unit.kilojoule_per_mole / omm_unit.nanometer**2,
        "kilojoule / (mole * nanometer ** 2)": omm_unit.kilojoule_per_mole / omm_unit.nanometer**2,
        "kj/mol/nm^2": omm_unit.kilojoule_per_mole / omm_unit.nanometer**2,
        # Temperature
        "kelvin": omm_unit.kelvin,
        "k": omm_unit.kelvin,
        # Pressure
        "bar": omm_unit.bar,
        "atmosphere": omm_unit.atmosphere,
        "atm": omm_unit.atmosphere,
        # Mass
        "dalton": omm_unit.dalton,
        "daltons": omm_unit.dalton,
        "amu": omm_unit.dalton,
        "gram": omm_unit.gram,
        "grams": omm_unit.gram,
        "kilogram": omm_unit.kilogram,
        "kilograms": omm_unit.kilogram,
        # Charge
        "elementary_charge": omm_unit.elementary_charge,
        "e": omm_unit.elementary_charge,
        # Concentration
        "molar": omm_unit.molar,
        "mole / liter": omm_unit.molar,
        "mol/l": omm_unit.molar,
        # Density
        "gram / milliliter": omm_unit.gram / omm_unit.milliliter,
        "g/ml": omm_unit.gram / omm_unit.milliliter,
        # Volume
        "nanometer ** 3": omm_unit.nanometer**3,
        "nm^3": omm_unit.nanometer**3,
        "angstrom ** 3": omm_unit.angstrom**3,
        "liter": omm_unit.liter,
        "liters": omm_unit.liter,
        "milliliter": omm_unit.milliliter,
        "milliliters": omm_unit.milliliter,
        # Dimensionless
        "dimensionless": omm_unit.dimensionless,
        "": omm_unit.dimensionless,
    }

    if unit_str in unit_map:
        return unit_map[unit_str]

    # Try parsing as a compound unit
    # This handles cases like "nanometer ** 3" that might have different spacing
    normalized = unit_str.replace(" ", "").replace("**", "^")

    # Check if normalized form matches
    for key, value in unit_map.items():
        if key.replace(" ", "").replace("**", "^") == normalized:
            return value

    raise ValueError(
        f"Cannot convert unit '{unit_str}' to OpenMM. "
        "Please add the unit mapping to polyzymd.utils.units._get_openmm_unit()"
    )


def _get_openff_unit_string(omm_unit_str: str) -> str:
    """Map an OpenMM unit string to an OpenFF/Pint unit string.

    Args:
        omm_unit_str: Unit string from OpenMM (e.g., "nanometer", "kilojoule/mole").

    Returns:
        The corresponding OpenFF unit string.

    Raises:
        ValueError: If the unit string is not recognized.
    """
    # Normalize the unit string
    omm_unit_str = omm_unit_str.strip()

    # Common unit mappings (OpenMM -> OpenFF)
    # OpenMM often uses different string representations
    unit_map = {
        # Length
        "nm": "nanometer",
        "nanometer": "nanometer",
        "angstrom": "angstrom",
        "A": "angstrom",
        "pm": "picometer",
        "picometer": "picometer",
        # Time
        "ps": "picosecond",
        "picosecond": "picosecond",
        "fs": "femtosecond",
        "femtosecond": "femtosecond",
        "ns": "nanosecond",
        "nanosecond": "nanosecond",
        # Energy
        "kJ/mol": "kilojoule/mole",
        "kilojoule/mole": "kilojoule/mole",
        "kcal/mol": "kilocalorie/mole",
        "kilocalorie/mole": "kilocalorie/mole",
        # Temperature
        "K": "kelvin",
        "kelvin": "kelvin",
        # Pressure
        "bar": "bar",
        "atm": "atmosphere",
        "atmosphere": "atmosphere",
        # Mass
        "Da": "dalton",
        "dalton": "dalton",
        "amu": "dalton",
        "g": "gram",
        "gram": "gram",
        "kg": "kilogram",
        "kilogram": "kilogram",
        # Charge
        "e": "elementary_charge",
        "elementary charge": "elementary_charge",
        # Volume
        "nm**3": "nanometer**3",
        "nm^3": "nanometer**3",
        "L": "liter",
        "liter": "liter",
        "mL": "milliliter",
        "milliliter": "milliliter",
        # Dimensionless
        "dimensionless": "dimensionless",
        "": "dimensionless",
    }

    if omm_unit_str in unit_map:
        return unit_map[omm_unit_str]

    # If not found, try returning as-is (Pint is often flexible)
    # This allows compound units to pass through
    return omm_unit_str


def ensure_openmm_units(
    value: Union[float, "OpenFFQuantity", "OpenMMQuantity"],
    default_unit: str = "nanometer",
) -> "OpenMMQuantity":
    """Ensure a value has OpenMM units, converting if necessary.

    This helper function accepts values in various forms:
    - Plain floats (assumes default_unit)
    - OpenFF Quantities (converts to OpenMM)
    - OpenMM Quantities (returns as-is)

    Args:
        value: A value that may or may not have units.
        default_unit: The unit to assume if value is a plain float.

    Returns:
        An OpenMM Quantity with units.

    Example:
        >>> # All of these return an OpenMM Quantity in nanometers:
        >>> ensure_openmm_units(1.5)  # float -> 1.5 nm
        >>> ensure_openmm_units(Quantity(1.5, "nm"))  # OpenFF -> OpenMM
        >>> ensure_openmm_units(1.5 * openmm.unit.nm)  # OpenMM -> OpenMM
    """
    import openmm.unit as omm_unit

    # Check if it's already an OpenMM Quantity
    if hasattr(value, "unit") and hasattr(value, "value_in_unit"):
        # This is an OpenMM Quantity - return as-is
        return value

    # Check if it's an OpenFF Quantity (has magnitude and units attributes)
    if hasattr(value, "magnitude") and hasattr(value, "units"):
        return openff_to_openmm(value)

    # It's a plain number - apply default unit
    default_omm_unit = _get_openmm_unit(default_unit)
    return float(value) * default_omm_unit


def ensure_openff_units(
    value: Union[float, "OpenFFQuantity", "OpenMMQuantity"],
    default_unit: str = "nanometer",
) -> "OpenFFQuantity":
    """Ensure a value has OpenFF units, converting if necessary.

    This helper function accepts values in various forms:
    - Plain floats (assumes default_unit)
    - OpenMM Quantities (converts to OpenFF)
    - OpenFF Quantities (returns as-is)

    Args:
        value: A value that may or may not have units.
        default_unit: The unit to assume if value is a plain float.

    Returns:
        An OpenFF Quantity with units.

    Example:
        >>> # All of these return an OpenFF Quantity in nanometers:
        >>> ensure_openff_units(1.5)  # float -> 1.5 nm
        >>> ensure_openff_units(1.5 * openmm.unit.nm)  # OpenMM -> OpenFF
        >>> ensure_openff_units(Quantity(1.5, "nm"))  # OpenFF -> OpenFF
    """
    from openff.units import Quantity as OFFQuantity

    # Check if it's an OpenFF Quantity (has magnitude and units attributes)
    if hasattr(value, "magnitude") and hasattr(value, "units"):
        return value

    # Check if it's an OpenMM Quantity
    if hasattr(value, "unit") and hasattr(value, "value_in_unit"):
        return openmm_to_openff(value)

    # It's a plain number - apply default unit
    return OFFQuantity(float(value), default_unit)
