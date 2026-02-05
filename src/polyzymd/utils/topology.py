"""
OpenFF Topology utilities.

This module provides functions for working with OpenFF Topology objects,
including loading from SDF files and extracting molecules.

The topology_from_sdf and related functions are adapted from the Polymerist
package by Timotej Bernat, used under the MIT License.

Original source: https://github.com/timbernat/polymerist
Copyright (c) 2024 Timotej Bernat

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import logging
from ast import literal_eval
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Union

if TYPE_CHECKING:
    from openff.toolkit import Molecule, Topology

LOGGER = logging.getLogger(__name__)


def _asiterable(obj: Any) -> Iterator:
    """Convert object to iterable, handling singleton case gracefully.

    If obj is already iterable (but not a string), returns it as-is.
    Otherwise, wraps it in a list.
    """
    if isinstance(obj, (list, tuple)):
        return iter(obj)
    elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes)):
        return iter(obj)
    else:
        return iter([obj])


def _package_atom_metadata(offmol: "Molecule") -> "Molecule":
    """Collect atom metadata into a serializable property of the parent Molecule.

    Creates a copy of the molecule with metadata packaged for serialization.

    Args:
        offmol: OpenFF Molecule to package metadata from.

    Returns:
        Copy of the molecule with metadata packaged into properties dict.
    """
    from copy import deepcopy

    packaged_mol = deepcopy(offmol)
    packaged_mdat: Dict[int, Dict] = {
        atom.molecule_atom_index: dict(atom.metadata)
        for atom in packaged_mol.atoms
        if atom.metadata
    }

    if packaged_mdat:  # only assign if metadata is present
        packaged_mol.properties["metadata"] = packaged_mdat

    return packaged_mol


def _unpackage_atom_metadata(offmol: "Molecule") -> "Molecule":
    """Reassign atom metadata from a packaged metadata property dict.

    Creates a copy of the molecule with metadata unpacked to atoms.

    Args:
        offmol: OpenFF Molecule with packaged metadata in properties.

    Returns:
        Copy of the molecule with metadata reassigned to atoms.
    """
    from copy import deepcopy

    unpacked_mol = deepcopy(offmol)
    packaged_mdat = unpacked_mol.properties.pop("metadata", None)

    if not packaged_mdat:  # handles both None and empty dict
        return unpacked_mol

    if isinstance(packaged_mdat, str):
        packaged_mdat = literal_eval(packaged_mdat)  # de-stringify if necessary

    for atom_id, metadata in packaged_mdat.items():
        # atom_id may be string after deserialization
        atom_idx = int(atom_id)
        unpacked_mol.atoms[atom_idx].metadata.update(metadata)

    return unpacked_mol


def topology_from_sdf(
    sdf_path: Union[str, Path],
    allow_undefined_stereo: bool = True,
    **kwargs: Any,
) -> "Topology":
    """Load an OpenFF Topology from an SDF file.

    Loads all molecules from the SDF file and creates a Topology
    containing them. Multi-molecule SDF files are supported.
    Atom metadata and partial charges stored in the SDF are preserved.

    This function is adapted from the Polymerist package by Timotej Bernat.

    Args:
        sdf_path: Path to the SDF file.
        allow_undefined_stereo: Whether to allow undefined stereochemistry.
        **kwargs: Additional arguments passed to Molecule.from_file().

    Returns:
        OpenFF Topology containing all molecules from the file.

    Raises:
        FileNotFoundError: If the SDF file doesn't exist.
        ValueError: If the file is not an SDF file.

    Example:
        >>> from polyzymd.utils.topology import topology_from_sdf
        >>> topology = topology_from_sdf("polymer.sdf")
        >>> print(f"Loaded {topology.n_molecules} molecules")
    """
    from openff.toolkit import Molecule, Topology

    sdf_path = Path(sdf_path)
    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    if sdf_path.suffix.lower() != ".sdf":
        raise ValueError(f"Expected .sdf file, got: {sdf_path.suffix}")

    LOGGER.debug(f"Loading serialized SDF Topology from {sdf_path}")

    # Load molecules and unpackage any serialized metadata
    molecules = Molecule.from_file(
        str(sdf_path),
        allow_undefined_stereo=allow_undefined_stereo,
        **kwargs,
    )

    # Build topology with metadata unpacked
    return Topology.from_molecules(_unpackage_atom_metadata(mol) for mol in _asiterable(molecules))


def topology_to_sdf(
    sdf_path: Union[str, Path],
    topology: "Topology",
) -> None:
    """Save an OpenFF Topology to an SDF file.

    Preserves all atom metadata by packaging it into molecule properties
    before serialization.

    This function is adapted from the Polymerist package by Timotej Bernat.

    Args:
        sdf_path: Path to save the SDF file.
        topology: OpenFF Topology to save.

    Raises:
        ValueError: If the path doesn't have .sdf extension.

    Example:
        >>> from polyzymd.utils.topology import topology_to_sdf
        >>> topology_to_sdf("output.sdf", my_topology)
    """
    sdf_path = Path(sdf_path)
    if sdf_path.suffix.lower() != ".sdf":
        raise ValueError(f"Expected .sdf extension, got: {sdf_path.suffix}")

    with sdf_path.open("w") as sdf_file:
        for mol in topology.molecules:
            # Package metadata for serialization without modifying original
            serial_mol = _package_atom_metadata(mol)
            serial_mol.to_file(sdf_file, file_format="SDF")

    LOGGER.debug(f"Successfully serialized SDF Topology to {sdf_path}")


def get_largest_offmol(topology: "Topology") -> "Molecule":
    """Extract the largest molecule from an OpenFF Topology.

    This is useful when loading SDF files that may contain multiple
    molecules (e.g., a polymer chain plus counterions). The largest
    molecule by atom count is typically the molecule of interest.

    This function is adapted from the Polymerist package by Timotej Bernat.

    Args:
        topology: OpenFF Topology containing one or more molecules.

    Returns:
        The molecule with the most atoms.

    Raises:
        ValueError: If the topology contains no molecules.

    Example:
        >>> from polyzymd.utils.topology import topology_from_sdf, get_largest_offmol
        >>> topology = topology_from_sdf("polymer_with_ions.sdf")
        >>> polymer = get_largest_offmol(topology)
        >>> print(f"Polymer has {polymer.n_atoms} atoms")
    """
    molecules: List["Molecule"] = list(topology.molecules)

    if not molecules:
        raise ValueError("Topology contains no molecules")

    # Find molecule with most atoms
    return max(molecules, key=lambda mol: mol.n_atoms)
