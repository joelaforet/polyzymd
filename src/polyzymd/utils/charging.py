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

    The default model is trained on AM1-BCC charges but can be customized.

    Attributes:
        model_name: Name/path of the NAGL model to use.

    Example:
        >>> charger = NAGLCharger()
        >>> charged_mol = charger.charge_molecule(molecule)
    """

    method_name: str = "nagl"

    def __init__(self, model_name: str = "openff-gnn-am1bcc-0.1.0-rc.2.pt") -> None:
        """Initialize the NAGL charger.

        Args:
            model_name: Name of the NAGL model. Can be:
                - A registered model name (e.g., "openff-gnn-am1bcc-0.1.0-rc.2.pt")
                - A path to a custom model file
        """
        self.model_name = model_name
        self._model = None

    def _load_model(self):
        """Lazily load the NAGL model.

        This defers the import and model loading until first use,
        allowing the module to be imported even if NAGL is not installed.
        """
        if self._model is None:
            try:
                from openff.nagl import GNNModel

                self._model = GNNModel.load(self.model_name)
            except ImportError as e:
                raise ImportError(
                    "NAGL is not installed. Install with: "
                    "mamba install -c conda-forge openff-nagl openff-nagl-models"
                ) from e
            except Exception as e:
                raise RuntimeError(f"Failed to load NAGL model '{self.model_name}': {e}") from e

        return self._model

    def charge_molecule(self, molecule: "Molecule") -> "Molecule":
        """Assign NAGL partial charges to a molecule.

        Args:
            molecule: OpenFF Molecule to charge.

        Returns:
            The molecule with NAGL charges assigned.

        Raises:
            ImportError: If NAGL is not installed.
            RuntimeError: If charging fails.
        """
        model = self._load_model()

        # NAGL doesn't require conformers but we generate one anyway
        # for consistency with other methods
        self._ensure_conformer(molecule)

        try:
            # NAGL assigns charges in-place
            model.compute_property(molecule, as_numpy=True)
        except Exception as e:
            raise RuntimeError(f"NAGL charge assignment failed: {e}") from e

        return molecule


class EspalomaCharger(MoleculeCharger):
    """Charger using Espaloma neural network.

    Espaloma (Extensible Surrogate Potential Optimized by Message-Passing
    Architectures) uses graph neural networks trained on quantum chemical
    data to assign partial charges and force field parameters.

    Attributes:
        model_name: Name/URL of the Espaloma model to use.

    Example:
        >>> charger = EspalomaCharger()
        >>> charged_mol = charger.charge_molecule(molecule)
    """

    method_name: str = "espaloma"

    def __init__(self, model_name: str = "espaloma-0.3.2") -> None:
        """Initialize the Espaloma charger.

        Args:
            model_name: Name of the Espaloma model. Can be:
                - A registered model version (e.g., "espaloma-0.3.2")
                - A URL to a model file
                - A local path to a model file
        """
        self.model_name = model_name
        self._model = None

    def _load_model(self):
        """Lazily load the Espaloma model.

        This defers the import and model loading until first use,
        allowing the module to be imported even if Espaloma is not installed.
        """
        if self._model is None:
            try:
                import espaloma as esp

                # Load model - espaloma handles URL/name resolution
                self._model = esp.get_model(self.model_name)
            except ImportError as e:
                raise ImportError(
                    "Espaloma is not installed. Install with: mamba install -c conda-forge espaloma"
                ) from e
            except Exception as e:
                raise RuntimeError(f"Failed to load Espaloma model '{self.model_name}': {e}") from e

        return self._model

    def charge_molecule(self, molecule: "Molecule") -> "Molecule":
        """Assign Espaloma partial charges to a molecule.

        Args:
            molecule: OpenFF Molecule to charge.

        Returns:
            The molecule with Espaloma charges assigned.

        Raises:
            ImportError: If Espaloma is not installed.
            RuntimeError: If charging fails.
        """
        try:
            import espaloma as esp
            from openff.units import Quantity
        except ImportError as e:
            raise ImportError(
                "Espaloma is not installed. Install with: mamba install -c conda-forge espaloma"
            ) from e

        model = self._load_model()

        # Espaloma requires a conformer
        self._ensure_conformer(molecule)

        try:
            # Create Espaloma graph from molecule
            mol_graph = esp.Graph(molecule)

            # Run model inference
            model(mol_graph.heterograph)

            # Extract charges from the graph
            # Espaloma stores charges in the 'n1' node type under 'q' key
            charges = mol_graph.heterograph.nodes["n1"].data["q"].flatten().numpy()

            # Assign charges to molecule
            molecule.partial_charges = Quantity(charges, "elementary_charge")

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
