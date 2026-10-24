"""The references a build record lists under "Methods to cite".

Each function takes the config and the build record (``engine``, the package
versions, ``provenance`` and ``components`` of ``build_provenance.json``) and
returns the references of one part of the build, or an empty list when the
build did not use it. :data:`CITATION_SOURCES` lists them in the order the record prints
them; :func:`polyzymd.simulation.artifact_integrity.methods_to_cite` joins
their results. To cite one more method, add a function and list it there.

Every DOI below was checked against Crossref (papers) or DataCite (Zenodo).
OpenFF parameters and software are cited as
https://openforcefield.org/science/how-to-cite/ asks: every reference it
lists for an item, the paper and the Zenodo DOI of the release.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable

FF14SB = (
    "Maier, J. A., Martinez, C., Kasavajhala, K., Wickstrom, L., Hauser, K. E. and "
    "Simmerling, C. (2015). ff14SB: improving the accuracy of protein side chain and "
    "backbone parameters from ff99SB. J. Chem. Theory Comput. 11(8), 3696-3713. "
    "doi:10.1021/acs.jctc.5b00255"
)
SAGE = (
    "Boothroyd, S., Behara, P. K., Madin, O. C., Hahn, D. F., Jang, H., Gapsys, V., "
    "Wagner, J. R., Horton, J. T., Dotson, D. L., Thompson, M. W., Maat, J., Gokey, T., "
    "Wang, L.-P., Cole, D. J., Gilson, M. K., Chodera, J. D., Bayly, C. I., Shirts, M. R. and "
    "Mobley, D. L. (2023). Development and benchmarking of Open Force Field 2.0.0: the Sage "
    "small molecule force field. J. Chem. Theory Comput. 19(11), 3251-3275. "
    "doi:10.1021/acs.jctc.3c00039"
)
PARSLEY = (
    "Qiu, Y., Smith, D. G. A., Boothroyd, S., Jang, H., Hahn, D. F., Wagner, J., Bannan, C. C., "
    "Gokey, T., Lim, V. T., Stern, C. D., Rizzi, A., Tjanaka, B., Tresadern, G., Lucas, X., "
    "Shirts, M. R., Gilson, M. K., Chodera, J. D., Bayly, C. I., Mobley, D. L. and Wang, L.-P. "
    "(2021). Development and benchmarking of Open Force Field v1.0.0, the Parsley "
    "small-molecule force field. J. Chem. Theory Comput. 17(10), 6262-6280. "
    "doi:10.1021/acs.jctc.1c00571"
)
SAGE_2_3 = (
    "Wang, L., Alibay, I., Behara, P. K., Boothroyd, S., Cavender, C. E., Horton, J. T., "
    "McIsaac, A. R., Mitchell, A., Morales, B., Thompson, M. W., Wagner, J. R., Westbrook, B. R., "
    "Bayly, C. I., Chodera, J. D., Cole, D. J., Eastwood, J. R. B., Shirts, M. R. and "
    "Mobley, D. L. (2026). Developing and benchmarking Sage 2.3.0 with the AshGC neural network "
    "charge model. J. Chem. Theory Comput. 22(9), 4507-4531. doi:10.1021/acs.jctc.6c00169"
)
TIP3P = (
    "Jorgensen, W. L., Chandrasekhar, J., Madura, J. D., Impey, R. W. and Klein, M. L. "
    "(1983). Comparison of simple potential functions for simulating liquid water. "
    "J. Chem. Phys. 79(2), 926-935. doi:10.1063/1.445869"
)
AM1BCC = (
    "Jakalian, A., Bush, B. L., Jack, D. B. and Bayly, C. I. (2000). Fast, efficient "
    "generation of high-quality atomic charges. AM1-BCC model: I. Method. J. Comput. Chem. "
    "21(2), 132-146. doi:10.1002/(SICI)1096-987X(20000130)21:2<132::AID-JCC5>3.0.CO;2-P",
    "Jakalian, A., Jack, D. B. and Bayly, C. I. (2002). Fast, efficient generation of "
    "high-quality atomic charges. AM1-BCC model: II. Parameterization and validation. "
    "J. Comput. Chem. 23(16), 1623-1641. doi:10.1002/jcc.10128",
)
ESPALOMA_CHARGE = (
    "Wang, Y., Pulido, I., Takaba, K., Kaminow, B., Scheen, J., Wang, L. and Chodera, J. D. "
    "(2024). EspalomaCharge: machine learning-enabled ultrafast partial charge assignment. "
    "J. Phys. Chem. A 128(20), 4160-4167. doi:10.1021/acs.jpca.4c01287"
)
JOUNG_CHEATHAM = (
    "Joung, I. S. and Cheatham, T. E. (2008). Determination of alkali and halide monovalent ion "
    "parameters for use in explicitly solvated biomolecular simulations. J. Phys. Chem. B "
    "112(30), 9020-9041. doi:10.1021/jp8001614"
)
OPENMM = (
    "Eastman, P., Galvelis, R., Peláez, R. P., Abreu, C. R. A., Farr, S. E., Gallicchio, E., "
    "Gorenko, A., Henry, M. M., Hu, F., Huang, J., Krämer, A., Michel, J., Mitchell, J. A., "
    "Pande, V. S., Rodrigues, J. P. G. L. M., Rodriguez-Guerra, J., Simmonett, A. C., Singh, S., "
    "Swails, J., Turner, P., Wang, Y., Zhang, I., Chodera, J. D., De Fabritiis, G. and "
    "Markland, T. E. (2024). OpenMM 8: molecular dynamics simulation with machine learning "
    "potentials. J. Phys. Chem. B 128(1), 109-116. doi:10.1021/acs.jpcb.3c06662"
)
GROMACS = (
    "Abraham, M. J., Murtola, T., Schulz, R., Páll, S., Smith, J. C., Hess, B. and Lindahl, E. "
    "(2015). GROMACS: high performance molecular simulations through multi-level parallelism "
    "from laptops to supercomputers. SoftwareX 1-2, 19-25. doi:10.1016/j.softx.2015.06.001"
)
OPENFF_TOOLKIT = (
    "Mobley, D. L., Bannan, C. C., Rizzi, A., Bayly, C. I., Chodera, J. D., Lim, V. T., "
    "Lim, N. M., Beauchamp, K. A., Slochower, D. R., Shirts, M. R., Gilson, M. K. and "
    "Eastman, P. K. (2018). Escaping atom types in force fields using direct chemical "
    "perception. J. Chem. Theory Comput. 14(11), 6076-6092. doi:10.1021/acs.jctc.8b00640"
)
OPENFF_INITIATIVE = (
    "Wang, L., Behara, P. K., Thompson, M. W., Gokey, T., Wang, Y., Wagner, J. R., Cole, D. J., "
    "Gilson, M. K., Shirts, M. R. and Mobley, D. L. (2024). The Open Force Field Initiative: "
    "open software and open science for molecular modeling. J. Phys. Chem. B 128(29), "
    "7043-7067. doi:10.1021/acs.jpcb.4c01558"
)
#: Zenodo DOIs OpenFF lists for its releases (each resolves to the release list).
OPENFF_TOOLKIT_ZENODO = "doi:10.5281/zenodo.597754"
OPENFF_INTERCHANGE_ZENODO = "doi:10.5281/zenodo.8147764"
OPENFF_FORCEFIELDS_ZENODO = "doi:10.5281/zenodo.3473553"
NAGL_REPOSITORY = "https://github.com/openforcefield/openff-nagl-models"


def _solutes(config: Any) -> list[Any]:
    """The solutes the build placed (a polymer entry can be disabled)."""
    return [s for s in getattr(config, "solutes", None) or [] if getattr(s, "enabled", True)]


def _co_solvents(config: Any) -> list[Any]:
    return list(getattr(getattr(config, "solvent", None), "co_solvents", None) or [])


def _force_field_files(config: Any) -> list[str]:
    """The force-field files the build loaded, in the order the config names them."""
    files: list[str] = []
    named = [getattr(s, "force_field", None) or [] for s in _solutes(config)]
    named.append(getattr(getattr(config, "solvent", None), "force_field", None) or [])
    for group in named:
        for name in group:
            if str(name) not in files:
                files.append(str(name))
    return files


def _openff_release(name: str) -> tuple[int, int] | None:
    """The (major, minor) release of an OpenFF Sage or Parsley file, by its name."""
    match = re.match(r"openff(?:_unconstrained)?-(\d+)\.(\d+)\.", Path(name).name.lower())
    return (int(match.group(1)), int(match.group(2))) if match else None


def _force_field_references(name: str) -> list[str]:
    """The references of one force-field file, by the file's name."""
    file = Path(name).name
    release = _openff_release(name)
    if "ff14sb" in file.lower():
        return [f"{file} (the OpenFF SMIRNOFF port of Amber ff14SB): {FF14SB}"]
    release_doi = f"{file} (the openff-forcefields release): {OPENFF_FORCEFIELDS_ZENODO}"
    if release is not None and release[0] == 2:
        cited = [f"{file} (a release of OpenFF Sage 2): {SAGE}"]
        if release == (2, 3):
            cited.append(f"{file} (OpenFF Sage 2.3.0): {SAGE_2_3}")
        return cited + [release_doi]
    if release is not None and release[0] == 1:
        return [f"{file} (a release of OpenFF Parsley 1): {PARSLEY}", release_doi]
    return [f"{file}: PolyzyMD has no reference for this file; cite its source."]


def cite_polyzymd(config: Any, record: dict[str, Any]) -> list[str]:
    """PolyzyMD itself, at the installed version."""
    from polyzymd.citation import DOI, citation_line

    return [citation_line() + (f" doi:{DOI}" if DOI else "")]


def cite_engine(config: Any, record: dict[str, Any]) -> list[str]:
    """The MD engine the build is for: OpenMM or GROMACS."""
    engine = record.get("engine")
    if engine == "openmm":
        version = record.get("openmm_version")
        return [f"OpenMM{f' {version}' if version else ''}: {OPENMM}"]
    if engine == "gromacs":
        return [f"GROMACS: {GROMACS}"]
    return []


def cite_openff_software(config: Any, record: dict[str, Any]) -> list[str]:
    """The OpenFF Toolkit, OpenFF Interchange and the Open Force Field Initiative.

    The Toolkit and Interchange parameterize every build. OpenFF asks for the
    Toolkit's paper and release DOI, Interchange's release DOI (it has no
    paper), and the Initiative's paper.
    """

    def release(key: str, doi: str) -> str:
        version = record.get(key)
        return f"{f'version {version}, ' if version else ''}{doi}"

    return [
        f"OpenFF Toolkit: {OPENFF_TOOLKIT}",
        f"OpenFF Toolkit release: {release('openff_toolkit_version', OPENFF_TOOLKIT_ZENODO)}",
        "OpenFF Interchange release: "
        f"{release('openff_interchange_version', OPENFF_INTERCHANGE_ZENODO)}",
        f"Open Force Field Initiative: {OPENFF_INITIATIVE}",
    ]


def cite_force_fields(config: Any, record: dict[str, Any]) -> list[str]:
    """Each force-field file the build loaded."""
    return [ref for name in _force_field_files(config) for ref in _force_field_references(name)]


def cite_water_model(config: Any, record: dict[str, Any]) -> list[str]:
    """The water model (PolyzyMD builds TIP3P only)."""
    model = getattr(getattr(getattr(config, "solvent", None), "water", None), "model", None)
    if getattr(model, "value", model) == "tip3p":
        return [f"TIP3P water: {TIP3P}"]
    return []


def cite_charge_models(config: Any, record: dict[str, Any]) -> list[str]:
    """The charge models of the molecules the build charged (NAGL, AM1-BCC, espaloma).

    Library charges belong to the force field cited above, and preset charges
    to the SDF files they came with.
    """
    specs = [getattr(s, "charges", None) for s in _solutes(config) + _co_solvents(config)]
    specs = [spec for spec in specs if spec is not None]
    methods = {spec.method for spec in specs}
    cited: list[str] = []
    for model in sorted({spec.model for spec in specs if spec.method == "nagl" and spec.model}):
        if Path(model).name.lower().startswith("openff-gnn-am1bcc-1."):
            cited.append(
                f"AshGC (OpenFF NAGL model {model}), a graph neural network for AM1-BCC "
                f"charges: {SAGE_2_3}"
            )
        else:
            cited.append(
                f"OpenFF NAGL charge model {model}, a graph neural network for AM1-BCC charges "
                f"({NAGL_REPOSITORY}); it has no journal article, so cite the model release and "
                "AM1-BCC."
            )
    shipped = sorted(
        {
            str(cosolvent.name)
            for cosolvent in _co_solvents(config)
            if getattr(getattr(cosolvent, "charges", None), "method", None) == "preset"
        }
    )
    if shipped:
        # scripts/generate_solvent_sdfs.py made the shipped charges with the
        # OpenFF Toolkit's AM1-BCC (docs: reference/configuration.md).
        cited.append(
            f"Charges of {', '.join(shipped)}: shipped with PolyzyMD, computed with AM1-BCC "
            "(cited below)."
        )
    if shipped or "nagl" in methods or "am1bcc" in methods:
        cited.extend(AM1BCC)
    if "espaloma" in methods:
        cited.append(ESPALOMA_CHARGE)
    return cited


def cite_ion_parameters(config: Any, record: dict[str, Any]) -> list[str]:
    """Joung and Cheatham, when the build added Na+ or Cl- with an OpenFF Sage or Parsley file.

    The Na+ and Cl- Lennard-Jones parameters of ``openff-1.*`` and ``openff-2.*``
    are the TIP3P set of Joung and Cheatham (2008), Table 5; the
    openff-forcefields package keeps that table as ``data/jc.csv``.
    """
    files = sorted(
        {
            Path(name).name
            for component in record.get("components") or []
            if component.get("type") == "ions"
            for domain in component.get("domains") or []
            for name in domain.get("force_field") or []
            if _openff_release(name) is not None and _openff_release(name)[0] in (1, 2)
        }
    )
    if not files:
        return []
    return [f"Na+ and Cl- parameters of {', '.join(files)} (TIP3P set): {JOUNG_CHEATHAM}"]


def cite_packmol(config: Any, record: dict[str, Any]) -> list[str]:
    """Packmol, when the build packed solvent or polymer chains with it (its seeds say so)."""
    from polyzymd.utils.packmol import PACKMOL

    provenance = record.get("provenance") or {}
    packed = ("solvent_packmol_seed", "polymer_packmol_seed")
    return [PACKMOL] if any(key in provenance for key in packed) else []


def cite_ion_rule(config: Any, record: dict[str, Any]) -> list[str]:
    """Machado and Pantano, when ``ionic_strength`` set the ions."""
    from polyzymd.builders.solvent import MACHADO_PANTANO

    ions = getattr(getattr(config, "solvent", None), "ions", None)
    return [MACHADO_PANTANO] if getattr(ions, "ionic_strength", None) is not None else []


def _references() -> list[str]:
    """Every full reference this module defines, longest first."""
    found: list[str] = []
    for value in globals().values():
        items = value if isinstance(value, tuple) else (value,)
        found += [item for item in items if isinstance(item, str) and " doi:" in item]
    return sorted(set(found), key=len, reverse=True)


def merge_repeated_references(cited: list[str]) -> list[str]:
    """Drop repeated references, keeping the first line that cites each one.

    A line is either a reference or ``"<what>: <reference>"``. When a later
    line cites the same reference for something else, its ``<what>`` joins the
    first line, so the reference is printed once and still says what it is for.
    """
    references = _references()
    merged: list[str] = []
    first: dict[str, int] = {}
    for line in cited:
        reference = next(
            (ref for ref in references if line == ref or line.endswith(": " + ref)), line
        )
        if reference not in first:
            first[reference] = len(merged)
            merged.append(line)
            continue
        what = line[: -len(reference) - 2] if line != reference else ""
        if what:
            index = first[reference]
            head = merged[index][: -len(reference) - 2] if merged[index] != reference else ""
            merged[index] = f"{head} and {what}: {reference}" if head else line
    return merged


#: The parts of a build that have references, in the order the record lists them.
CITATION_SOURCES: tuple[Callable[[Any, dict[str, Any]], list[str]], ...] = (
    cite_polyzymd,
    cite_engine,
    cite_openff_software,
    cite_force_fields,
    cite_water_model,
    cite_ion_parameters,
    cite_charge_models,
    cite_packmol,
    cite_ion_rule,
)
