"""The references a build record lists under "Methods to cite".

Each reference is one :class:`Reference` in :data:`REFERENCES`, in print
order: the reference as printed (authors, year, title, journal, DOI), when the
build used what it describes (a predicate of the :class:`BuildFacts` read once
from the config and the build record) and what it is cited for. To cite one
more method, add one entry, e.g. ``Reference(PAPER, "mBuild", lambda f, _:
"polymer_builder" in f.provenance)``; :func:`cite` lists it.

Every DOI below was checked against Crossref (papers) or DataCite (Zenodo).
OpenFF parameters and software are cited as
https://openforcefield.org/science/how-to-cite/ asks: every reference it
lists for an item, the paper and the Zenodo DOI of the release.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from polyzymd import citation
from polyzymd.builders.polymer_generator import MBUILD
from polyzymd.builders.solvent import MACHADO_PANTANO
from polyzymd.utils.packmol import PACKMOL

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


def _release(name: Any) -> tuple[int, int]:
    """(major, minor) of an OpenFF Parsley (1) or Sage (2) file, (-1, 0) for ff14SB, else (0, 0)."""
    file = Path(name).name.lower()
    match = re.match(r"openff(?:_unconstrained)?-([12])\.(\d+)\.", file)
    if "ff14sb" in file:
        return (-1, 0)
    return (int(match[1]), int(match[2])) if match else (0, 0)


def _joined(names: Any) -> str:
    return ", ".join(sorted(set(names)))


@dataclass(frozen=True)
class BuildFacts:
    """What a build used, read once from the config and the record for :data:`REFERENCES`."""

    record: dict[str, Any]
    provenance: dict[str, Any]
    force_fields: tuple[Path, ...]  # each file once, in the order the config names them
    water: Any
    charge_methods: frozenset[str]  # with "am1bcc" for NAGL and shipped charges, AM1-BCC too
    nagl_models: tuple[str, ...]
    shipped_co_solvents: str  # co-solvents with charges shipped with PolyzyMD, joined
    ion_files: str  # the OpenFF files that gave the ions their parameters, joined
    ionic_strength: Any
    packmol: bool  # whether the build packed solvent or polymer chains (its seeds say so)

    @classmethod
    def of(cls, config: Any, record: dict[str, Any]) -> BuildFacts:
        solvent = getattr(config, "solvent", None)
        # A polymer entry can be disabled; it is then not built.
        solutes = [s for s in getattr(config, "solutes", None) or [] if getattr(s, "enabled", True)]
        co_solvents = list(getattr(solvent, "co_solvents", None) or [])
        names = [n for s in solutes for n in getattr(s, "force_field", None) or []]
        names += getattr(solvent, "force_field", None) or []
        specs = [s.charges for s in solutes + co_solvents if getattr(s, "charges", None)]
        model = getattr(getattr(solvent, "water", None), "model", None)
        ions = [
            name
            for component in record.get("components") or []
            if component.get("type") == "ions"
            for domain in component.get("domains") or []
            for name in domain.get("force_field") or []
        ]
        shipped = _joined(
            str(c.name)
            for c in co_solvents
            if getattr(getattr(c, "charges", None), "method", None) == "preset"
        )
        methods = {spec.method for spec in specs}
        provenance = record.get("provenance") or {}
        return cls(
            record=record,
            provenance=provenance,
            force_fields=tuple(Path(n) for n in dict.fromkeys(str(n) for n in names)),
            water=getattr(model, "value", model),
            charge_methods=frozenset(
                methods | ({"am1bcc"} if shipped or "nagl" in methods else set())
            ),
            nagl_models=tuple(sorted({s.model for s in specs if s.method == "nagl" and s.model})),
            shipped_co_solvents=shipped,
            ion_files=_joined(Path(name).name for name in ions if _release(name)[0] > 0),
            ionic_strength=getattr(getattr(solvent, "ions", None), "ionic_strength", None),
            packmol=bool({"solvent_packmol_seed", "polymer_packmol_seed"} & provenance.keys()),
        )

    def release(self, package: str, doi: str) -> str:
        """``doi`` of an OpenFF package release, after the version the record names."""
        version = self.record.get(f"openff_{package}_version")
        return f"{f'version {version}, ' if version else ''}{doi}"


@dataclass(frozen=True)
class Reference:
    """One line of "Methods to cite": a reference, what it is cited for, and when.

    ``label`` is a format string of ``f`` (the facts) and ``x`` (the item), or
    a function ``(facts, item)``; ``text`` may be such a function for the few
    lines that carry a version or a model name. ``when(facts, item)`` says
    whether the build used it. ``each`` names a tuple of :class:`BuildFacts`
    to check it for, item by item; ``item`` is None without it.
    """

    text: str | Callable[[BuildFacts, Any], str]
    label: str | Callable[[BuildFacts, Any], str] = ""
    when: Callable[[BuildFacts, Any], Any] = lambda facts, item: True
    each: str = ""

    def line(self, facts: BuildFacts, item: Any) -> str:
        text = self.text if isinstance(self.text, str) else self.text(facts, item)
        label = self.label
        label = label.format(f=facts, x=item) if isinstance(label, str) else label(facts, item)
        return f"{label}: {text}" if label else text


def _file(text: str, what: str, when: Callable[[Path], Any]) -> Reference:
    """A reference cited for each force-field file ``when`` matches, by the file's name."""
    return Reference(text, "{x.name}" + what, lambda f, x: when(x), each="force_fields")


def _ashgc(facts: BuildFacts, model: str) -> bool:
    return Path(model).name.lower().startswith("openff-gnn-am1bcc-1.")


NO_REFERENCE = "PolyzyMD has no reference for this file; cite its source."
#: Every reference "Methods to cite" can list, in the order it lists them.
REFERENCES: tuple[Reference, ...] = (
    Reference(lambda f, _: " doi:".join(filter(None, (citation.citation_line(), citation.DOI)))),
    Reference(
        OPENMM,
        lambda f, _: " ".join(filter(None, ("OpenMM", f.record.get("openmm_version")))),
        lambda f, _: f.record.get("engine") == "openmm",
    ),
    Reference(GROMACS, "GROMACS", lambda f, _: f.record.get("engine") == "gromacs"),
    # OpenFF asks for the Toolkit's paper and release DOI, Interchange's release
    # DOI (it has no paper), and the Initiative's paper; both run every build.
    Reference(OPENFF_TOOLKIT, "OpenFF Toolkit"),
    Reference(lambda f, _: f.release("toolkit", OPENFF_TOOLKIT_ZENODO), "OpenFF Toolkit release"),
    Reference(
        lambda f, _: f.release("interchange", OPENFF_INTERCHANGE_ZENODO),
        "OpenFF Interchange release",
    ),
    Reference(OPENFF_INITIATIVE, "Open Force Field Initiative"),
    # The Sage paper describes 2.0.0 and is cited for every Sage 2 release.
    _file(FF14SB, " (the OpenFF SMIRNOFF port of Amber ff14SB)", lambda x: _release(x)[0] < 0),
    _file(SAGE, " (a release of OpenFF Sage 2)", lambda x: _release(x)[0] == 2),
    _file(SAGE_2_3, " (OpenFF Sage 2.3.0)", lambda x: _release(x) == (2, 3)),
    _file(PARSLEY, " (a release of OpenFF Parsley 1)", lambda x: _release(x)[0] == 1),
    _file(
        OPENFF_FORCEFIELDS_ZENODO, " (the openff-forcefields release)", lambda x: _release(x)[0] > 0
    ),
    _file(NO_REFERENCE, "", lambda x: _release(x)[0] == 0),
    Reference(TIP3P, "TIP3P water", lambda f, _: f.water == "tip3p"),
    # The Na+ and Cl- Lennard-Jones parameters of openff-1.* and openff-2.* are
    # the TIP3P set of Joung and Cheatham (2008), Table 5 (data/jc.csv of
    # openff-forcefields).
    Reference(
        JOUNG_CHEATHAM,
        "Na+ and Cl- parameters of {f.ion_files} (TIP3P set)",
        lambda f, _: f.ion_files,
    ),
    # Charge models. Library charges belong to the force field cited above, and
    # preset charges to the SDF files they came with.
    Reference(
        SAGE_2_3,
        "AshGC (OpenFF NAGL model {x}), a graph neural network for AM1-BCC charges",
        _ashgc,
        each="nagl_models",
    ),
    Reference(
        lambda f, model: (
            f"OpenFF NAGL charge model {model}, a graph neural network for AM1-BCC charges "
            f"({NAGL_REPOSITORY}); it has no journal article, so cite the model release and "
            "AM1-BCC."
        ),
        when=lambda f, x: not _ashgc(f, x),
        each="nagl_models",
    ),
    # scripts/generate_solvent_sdfs.py made the shipped charges with the
    # OpenFF Toolkit's AM1-BCC (docs: reference/configuration.md).
    Reference(
        "shipped with PolyzyMD, computed with AM1-BCC (cited below).",
        "Charges of {f.shipped_co_solvents}",
        lambda f, _: f.shipped_co_solvents,
    ),
    *(Reference(paper, when=lambda f, _: "am1bcc" in f.charge_methods) for paper in AM1BCC),
    Reference(ESPALOMA_CHARGE, when=lambda f, _: "espaloma" in f.charge_methods),
    Reference(PACKMOL, when=lambda f, _: f.packmol),
    Reference(MACHADO_PANTANO, when=lambda f, _: f.ionic_strength is not None),
    # mBuild, when the build generated polymer chains with it.
    Reference(MBUILD, when=lambda f, _: "polymer_builder" in f.provenance),
)


def cite(config: Any, record: dict[str, Any]) -> list[str]:
    """The lines of "Methods to cite" for a build, each reference once."""
    facts = BuildFacts.of(config, record)
    cited: list[str] = []
    for each, group in itertools.groupby(REFERENCES, key=lambda reference: reference.each):
        entries = list(group)
        for item in getattr(facts, each) if each else (None,):
            cited += [ref.line(facts, item) for ref in entries if ref.when(facts, item)]
    return merge_repeated_references(cited)


def _references() -> list[str]:
    """Every full reference of :data:`REFERENCES`, longest first."""
    found = {r.text for r in REFERENCES if isinstance(r.text, str) and " doi:" in r.text}
    return sorted(found, key=len, reverse=True)


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
