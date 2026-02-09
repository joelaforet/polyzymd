"""
Molecular charging utilities with extensible charger classes.

This module provides a Strategy pattern implementation for assigning partial
charges to molecules using various ML-based and QM-based methods:

- NAGLCharger: NAGL neural network charges (fast, graph-based)
- EspalomaCharger: Espaloma neural network charges
- AM1BCCCharger: Semi-empirical AM1-BCC charges (OpenFF toolkit built-in)

The abstract MoleculeCharger base class defines a consistent interface that
allows easy addition of future charging methods.

Example:
    >>> from openff.toolkit import Molecule
    >>> from polyzymd.utils.charging import get_charger
    >>>
    >>> mol = Molecule.from_smiles("CCO")  # Ethanol
    >>> charger = get_charger("nagl")
    >>> charged_mol = charger.charge_molecule(mol)
    >>> print(charged_mol.partial_charges)

Adapted from the Polymerist package by Timotej Bernat, used under the MIT License.
Original source: https://github.com/timbernat/polymerist
Copyright (c) 2024 Timotej Bernat
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional, Type

if TYPE_CHECKING:
    from openff.toolkit import Molecule


# Default NAGL model - using rc.3 to match Polymerist
# TODO: Update to production model "openff-gnn-am1bcc-1.0.0.pt" when ready
#       (requires openff-nagl-models >= v2025.09.0)
DEFAULT_NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"

# Default Espaloma charge method - this is currently the only method supported
# by EspalomaChargeToolkitWrapper, but we keep it configurable for future flexibility
DEFAULT_ESPALOMA_METHOD = "espaloma-am1bcc"


class MoleculeCharger(ABC):
    """Abstract base class for molecular charging strategies.

    Subclasses must implement the charge_molecule() method to assign
    partial charges using their specific methodology.

    This follows the Strategy pattern, allowing different charging
    algorithms to be used interchangeably.

    Attributes:
        method_name: Human-readable name of the charging method.
    """

    method_name: str = "base"

    @abstractmethod
    def charge_molecule(self, molecule: "Molecule") -> "Molecule":
        """Assign partial charges to a molecule.

        Args:
            molecule: OpenFF Molecule to charge. The molecule should have
                at least one conformer; if not, one will be generated.

        Returns:
            The same molecule with partial_charges assigned. Note that this
            modifies the molecule in-place AND returns it for convenience.

        Raises:
            RuntimeError: If charging fails for any reason.
        """
        pass

    def _ensure_conformer(self, molecule: "Molecule") -> None:
        """Ensure the molecule has at least one conformer.

        Many charging methods require 3D coordinates. This method generates
        a conformer using the OpenFF toolkit if none exists.

        Args:
            molecule: Molecule to check/modify.
        """
        if molecule.n_conformers == 0:
            molecule.generate_conformers(n_conformers=1)


class NAGLCharger(MoleculeCharger):
    """Charger using NAGL (Neural Network Atomic Potential Library).

    NAGL uses graph neural networks to predict partial charges directly
    from molecular graphs, making it very fast and suitable for large-scale
    screening applications.

    This implementation uses the NAGLToolkitWrapper from OpenFF Toolkit,
    which handles model path resolution via the openff-nagl-models package.

    Attributes:
        model_name: Name of the NAGL model to use. Must be a model name
            registered in openff-nagl-models (e.g., "openff-gnn-am1bcc-0.1.0-rc.3.pt").

    Example:
        >>> charger = NAGLCharger()
        >>> charged_mol = charger.charge_molecule(molecule)
    """

    method_name: str = "nagl"

    def __init__(self, model_name: str = DEFAULT_NAGL_MODEL) -> None:
        """Initialize the NAGL charger.

        Args:
            model_name: Name of the NAGL model. Should be a model name
                registered in the openff-nagl-models package, such as:
                - "openff-gnn-am1bcc-0.1.0-rc.3.pt" (default, pre-production)
                - "openff-gnn-am1bcc-1.0.0.pt" (production, requires newer package)
        """
        self.model_name = model_name

    def charge_molecule(self, molecule: "Molecule") -> "Molecule":
        """Assign NAGL partial charges to a molecule.

        Uses NAGLToolkitWrapper from OpenFF Toolkit, which handles
        model path resolution automatically via openff-nagl-models.

        Args:
            molecule: OpenFF Molecule to charge.

        Returns:
            The molecule with NAGL charges assigned.

        Raises:
            ImportError: If NAGL or nagl-models is not installed.
            RuntimeError: If charging fails.
        """
        try:
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        except ImportError as e:
            raise ImportError(
                "NAGL is not installed. Install with: "
                "mamba install -c conda-forge openff-nagl openff-nagl-models"
            ) from e

        # NAGL doesn't strictly require conformers but we generate one anyway
        # for consistency with other methods
        self._ensure_conformer(molecule)

        try:
            molecule.assign_partial_charges(
                partial_charge_method=self.model_name,
                toolkit_registry=NAGLToolkitWrapper(),
            )
        except Exception as e:
            raise RuntimeError(f"NAGL charge assignment failed: {e}") from e

        return molecule


class EspalomaCharger(MoleculeCharger):
    """Charger using Espaloma neural network.

    Espaloma (Extensible Surrogate Potential Optimized by Message-Passing
    Architectures) uses graph neural networks trained on quantum chemical
    data to assign partial charges and force field parameters.

    This implementation uses the EspalomaChargeToolkitWrapper from the
    espaloma-charge package, which provides a clean interface to the
    Espaloma charging functionality.

    Attributes:
        charge_method: The charge method to use with Espaloma.
            Currently only "espaloma-am1bcc" is supported by the toolkit wrapper.

    Example:
        >>> charger = EspalomaCharger()
        >>> charged_mol = charger.charge_molecule(molecule)
    """

    method_name: str = "espaloma"

    def __init__(self, charge_method: str = DEFAULT_ESPALOMA_METHOD) -> None:
        """Initialize the Espaloma charger.

        Args:
            charge_method: The charge method to use. Currently only
                "espaloma-am1bcc" is supported by EspalomaChargeToolkitWrapper,
                but this parameter is kept for future flexibility.
        """
        self.charge_method = charge_method

    def charge_molecule(self, molecule: "Molecule") -> "Molecule":
        """Assign Espaloma partial charges to a molecule.

        Uses EspalomaChargeToolkitWrapper from the espaloma-charge package,
        which provides a clean interface matching the OpenFF Toolkit pattern.

        Args:
            molecule: OpenFF Molecule to charge.

        Returns:
            The molecule with Espaloma charges assigned.

        Raises:
            ImportError: If espaloma-charge is not installed.
            RuntimeError: If charging fails.
        """
        try:
            from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
        except ImportError as e:
            raise ImportError(
                "espaloma-charge is not installed. Install with: "
                "mamba install -c conda-forge espaloma-charge"
            ) from e

        # Espaloma requires a conformer
        self._ensure_conformer(molecule)

        try:
            molecule.assign_partial_charges(
                partial_charge_method=self.charge_method,
                toolkit_registry=EspalomaChargeToolkitWrapper(),
            )
        except Exception as e:
            raise RuntimeError(f"Espaloma charge assignment failed: {e}") from e

        return molecule


class AM1BCCCharger(MoleculeCharger):
    """Charger using AM1-BCC semi-empirical method.

    AM1-BCC (Austin Model 1 with Bond Charge Corrections) computes charges
    using the AM1 semi-empirical quantum method with empirical corrections.

    This is the default charging method in the OpenFF toolkit and uses
    either AmberTools (sqm) or the OpenEye toolkit as the backend.

    Example:
        >>> charger = AM1BCCCharger()
        >>> charged_mol = charger.charge_molecule(molecule)
    """

    method_name: str = "am1bcc"

    def __init__(self, toolkit: Optional[str] = None) -> None:
        """Initialize the AM1-BCC charger.

        Args:
            toolkit: Which toolkit to use for AM1-BCC calculation.
                - None (default): Use the first available toolkit
                - "ambertools": Use AmberTools (sqm)
                - "openeye": Use OpenEye toolkit
        """
        self.toolkit = toolkit

    def charge_molecule(self, molecule: "Molecule") -> "Molecule":
        """Assign AM1-BCC partial charges to a molecule.

        Args:
            molecule: OpenFF Molecule to charge.

        Returns:
            The molecule with AM1-BCC charges assigned.

        Raises:
            RuntimeError: If charging fails.
        """
        # AM1-BCC requires a conformer
        self._ensure_conformer(molecule)

        try:
            if self.toolkit == "openeye":
                from openff.toolkit.utils import OpenEyeToolkitWrapper

                molecule.assign_partial_charges("am1bcc", toolkit_registry=OpenEyeToolkitWrapper())
            elif self.toolkit == "ambertools":
                from openff.toolkit.utils import AmberToolsToolkitWrapper

                molecule.assign_partial_charges(
                    "am1bcc", toolkit_registry=AmberToolsToolkitWrapper()
                )
            else:
                # Use default toolkit priority
                molecule.assign_partial_charges("am1bcc")

        except Exception as e:
            raise RuntimeError(f"AM1-BCC charge assignment failed: {e}") from e

        return molecule


# Registry of available chargers
_CHARGER_REGISTRY: dict[str, Type[MoleculeCharger]] = {
    "nagl": NAGLCharger,
    "espaloma": EspalomaCharger,
    "am1bcc": AM1BCCCharger,
}


def get_charger(method: str, **kwargs) -> MoleculeCharger:
    """Factory function to get a charger by method name.

    This is the recommended way to obtain charger instances, as it
    provides a consistent interface and handles method name resolution.

    Args:
        method: Name of the charging method. Case-insensitive.
            Supported values: "nagl", "espaloma", "am1bcc"
        **kwargs: Additional arguments passed to the charger constructor.

    Returns:
        An instance of the appropriate MoleculeCharger subclass.

    Raises:
        ValueError: If the method is not recognized.

    Example:
        >>> charger = get_charger("nagl")
        >>> charger = get_charger("am1bcc", toolkit="openeye")
    """
    method_lower = method.lower().strip()

    if method_lower not in _CHARGER_REGISTRY:
        available = ", ".join(sorted(_CHARGER_REGISTRY.keys()))
        raise ValueError(f"Unknown charging method: '{method}'. Available methods: {available}")

    charger_class = _CHARGER_REGISTRY[method_lower]
    return charger_class(**kwargs)


def register_charger(name: str, charger_class: Type[MoleculeCharger]) -> None:
    """Register a custom charger class.

    This allows users to add their own charging methods that integrate
    with the get_charger() factory function.

    Args:
        name: Name to register the charger under (case-insensitive).
        charger_class: A subclass of MoleculeCharger.

    Raises:
        TypeError: If charger_class is not a subclass of MoleculeCharger.

    Example:
        >>> class MyCustomCharger(MoleculeCharger):
        ...     method_name = "custom"
        ...     def charge_molecule(self, molecule):
        ...         # Custom charging logic
        ...         return molecule
        >>> register_charger("custom", MyCustomCharger)
        >>> charger = get_charger("custom")
    """
    if not isinstance(charger_class, type) or not issubclass(charger_class, MoleculeCharger):
        raise TypeError(
            f"charger_class must be a subclass of MoleculeCharger, got {type(charger_class)}"
        )

    _CHARGER_REGISTRY[name.lower().strip()] = charger_class
