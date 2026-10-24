"""The references a build record lists under "Methods to cite".

Each function takes the config and the build's provenance and returns the
references of one part of the build, or an empty list when the build did not
use it. :data:`CITATION_SOURCES` lists them in the order the record prints
them; :func:`polyzymd.simulation.artifact_integrity.methods_to_cite` joins
their results. To cite one more method, add a function and list it there.

Every DOI below was checked against Crossref.
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


def _force_field_reference(name: str) -> str:
    """The reference of one force-field file, by the file's name."""
    stem = Path(name).name.lower()
    if "ff14sb" in stem:
        return f"{Path(name).name} (the OpenFF SMIRNOFF port of Amber ff14SB): {FF14SB}"
    match = re.match(r"openff(?:_unconstrained)?-(\d+)\.", stem)
    if match and match.group(1) == "2":
        return f"{Path(name).name} (a release of OpenFF Sage 2): {SAGE}"
    if match and match.group(1) == "1":
        return f"{Path(name).name} (a release of OpenFF Parsley 1): {PARSLEY}"
    return f"{Path(name).name}: PolyzyMD has no reference for this file; cite its source."


def cite_polyzymd(config: Any, provenance: dict[str, Any]) -> list[str]:
    """PolyzyMD itself, at the installed version."""
    from polyzymd.citation import DOI, citation_line

    return [citation_line() + (f" doi:{DOI}" if DOI else "")]


def cite_force_fields(config: Any, provenance: dict[str, Any]) -> list[str]:
    """Each force-field file the build loaded."""
    return [_force_field_reference(name) for name in _force_field_files(config)]


def cite_water_model(config: Any, provenance: dict[str, Any]) -> list[str]:
    """The water model (PolyzyMD builds TIP3P only)."""
    model = getattr(getattr(getattr(config, "solvent", None), "water", None), "model", None)
    if getattr(model, "value", model) == "tip3p":
        return [f"TIP3P water: {TIP3P}"]
    return []


def cite_charge_models(config: Any, provenance: dict[str, Any]) -> list[str]:
    """The charge models of the molecules the build charged (NAGL, AM1-BCC, espaloma).

    Library charges belong to the force field cited above, and preset charges
    to the SDF files they came with.
    """
    specs = [getattr(s, "charges", None) for s in _solutes(config) + _co_solvents(config)]
    specs = [spec for spec in specs if spec is not None]
    methods = {spec.method for spec in specs}
    cited: list[str] = []
    for model in sorted({spec.model for spec in specs if spec.method == "nagl" and spec.model}):
        cited.append(
            f"OpenFF NAGL charge model {model}, a graph neural network for AM1-BCC charges "
            f"({NAGL_REPOSITORY}); it has no journal article, so cite the model release and "
            "AM1-BCC."
        )
    if "nagl" in methods or "am1bcc" in methods:
        cited.extend(AM1BCC)
    if "espaloma" in methods:
        cited.append(ESPALOMA_CHARGE)
    return cited


def cite_packmol(config: Any, provenance: dict[str, Any]) -> list[str]:
    """Packmol, when the build packed solvent or polymer chains with it (its seeds say so)."""
    from polyzymd.utils.packmol import PACKMOL

    packed = ("solvent_packmol_seed", "polymer_packmol_seed")
    return [PACKMOL] if any(key in provenance for key in packed) else []


def cite_ion_rule(config: Any, provenance: dict[str, Any]) -> list[str]:
    """Machado and Pantano, when ``ionic_strength`` set the ions."""
    from polyzymd.builders.solvent import MACHADO_PANTANO

    ions = getattr(getattr(config, "solvent", None), "ions", None)
    return [MACHADO_PANTANO] if getattr(ions, "ionic_strength", None) is not None else []


#: The parts of a build that have references, in the order the record lists them.
CITATION_SOURCES: tuple[Callable[[Any, dict[str, Any]], list[str]], ...] = (
    cite_polyzymd,
    cite_force_fields,
    cite_water_model,
    cite_charge_models,
    cite_packmol,
    cite_ion_rule,
)
