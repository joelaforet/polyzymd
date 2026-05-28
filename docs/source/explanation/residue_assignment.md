# Understanding Residue Assignment in PolyzyMD

PolyzyMD uses chain IDs and residue numbers to preserve chemically meaningful
identity in generated topology files. The goal is not only to produce a valid
simulation system, but also to keep enough identity information for selection,
tracking, visualization, and summary analyses after the system has been built.

This page explains the design convention behind those identifiers. It is meant
for contributors and advanced users who need to understand why generated
PolyzyMD topologies are organized the way they are.

## The core principle: one repeat unit is one residue

PolyzyMD treats the smallest chemically meaningful unit as the residue-level
unit whenever possible:

> **One repeat unit is one residue.**

For proteins, this matches the usual biomolecular meaning of a residue: each
amino acid is one residue. For synthetic polymers, the analogous unit is the
monomer or repeat unit. For solvent, ions, and other small molecules, the
smallest meaningful unit is usually the complete molecule.

This convention makes residue identifiers useful for analysis rather than just
file-format bookkeeping.

| Component type | Chemically meaningful unit | Residue assignment |
|----------------|----------------------------|--------------------|
| Protein | Amino acid | Each amino acid remains one residue. |
| Substrate or ligand | Complete molecule, unless represented otherwise | The ligand is assigned as a distinct unit. |
| Polymer | Monomer or source-defined repeat unit | Each represented repeat unit can be selected and summarized. |
| Solvent and ions | Complete molecule or ion | Each molecule or ion can be distinguished. |
| Co-solvent | Complete molecule | Each molecule can be distinguished. |

The practical consequence is that a residue identifier should point to a unit a
scientist might reasonably select, track, count, or summarize.

## Why default topology output can lose information

Some topology-generation workflows preserve atom names and coordinates but do
not preserve useful residue-level identity for every molecule. In the most
problematic case, many molecules of the same type can appear to share a single
residue identity.

Conceptually, this loses information in several ways:

- Distinct solvent molecules or ions cannot be cleanly distinguished by residue
  identity.
- Per-molecule summaries, such as residence time or hydrogen-bond counts, become
  harder to define.
- Visualization tools cannot reliably refer to a specific molecule or polymer
  unit by a stable identifier.
- Analysis results become less interpretable because atom-level selections are
  no longer tied to chemically meaningful units.

PolyzyMD chooses to preserve this identity in the generated topology because it
is difficult or impossible to reconstruct unambiguously after the identifiers
have been discarded.

## Chain convention

PolyzyMD uses a consistent chain convention in generated topology files:

| Chain | Role |
|-------|------|
| A | Protein or enzyme. |
| B | Substrate or ligand. |
| C | Polymer or conjugated polymer component. |
| D and later | Solvent, ions, and remaining molecules. |

The convention is intentionally simple. It gives contributors, analysis
plugins, and visualization workflows a shared vocabulary for referring to major
system components without needing to infer component roles from atom names or
force-field metadata.

For large systems, solvent and other remaining molecules may span multiple chain
IDs because common topology formats place limits on residue numbering within a
chain. In that case, the role is still the same: chains D and later represent
solvent, ions, or other non-protein, non-ligand, non-polymer components.

## Polymer identity and residue numbering

Polymer systems need extra care because the chemically meaningful unit may come
from the source representation. PolyzyMD's residue assignment should preserve
distinct monomer or repeat-unit identity when that identity is present.

This does not mean contributors should assume every polymer chain always
restarts residue numbering. Whether numbering restarts depends on the source
representation and the generated topology semantics. The invariant is stronger
and more important than a particular numbering pattern:

- Distinct represented polymer units should remain distinguishable.
- Chain and residue identifiers together should identify the intended unit.
- Source-defined monomer identity should not be flattened into a generic polymer
  residue when that would make per-unit analysis ambiguous.

For example, two polymer strands in the same system should preserve enough
identity to distinguish units on one strand from units on another. A contributor
should not write analysis code that assumes residue number 1 always means “the
first monomer of every polymer chain” unless the topology source explicitly
defines that convention.

## What the convention enables

The residue and chain convention is not an analysis method by itself. It is a
foundation that makes downstream methods easier to write and interpret.

Because generated topologies preserve chemically meaningful unit identity,
analyses can ask questions such as:

- Which solvent molecules remain near an active site across frames?
- Which polymer repeat units contact the protein most often?
- Which ligand or substrate atoms are involved in recurring interactions with
  chain A?
- How do per-monomer properties vary along chain C?
- Which waters, ions, or co-solvent molecules bridge between protein and ligand?

Visualization workflows benefit for the same reason: interesting molecules or
monomers identified by analysis can be mapped back to chain and residue
identifiers in the generated topology.

## Generated topology behavior

PolyzyMD-generated topology files expose the assignment convention through
standard chain and residue identifiers. Contributors should describe and rely on
that public behavior rather than private implementation details.

At a high level, generated topologies aim to preserve these properties:

- Protein residue identities from the input structure remain meaningful.
- The substrate or ligand is separated from the protein by chain identity.
- Polymer units on chain C remain distinguishable at residue granularity.
- Solvent, ions, and remaining molecules on chains D and later remain
  distinguishable as chemically meaningful units.
- The combined chain/residue identity is the stable handle for a unit when a
  residue number alone might be ambiguous.

This design prioritizes analysis capability over the smallest possible topology
metadata representation. A generated file may carry more residue identity
information than a minimal simulation input, but that information is valuable
for reproducible interpretation.

## Conceptual checklist for reviewing assignments

When reviewing generated topology behavior or changing builder code that affects
identifiers, use this checklist conceptually:

- Can a user distinguish the protein, ligand, polymer, and solvent components by
  chain identity?
- Does one residue correspond to one chemically meaningful repeat unit for the
  component being represented?
- Are solvent molecules, ions, and co-solvents distinguishable as individual
  units when analysis needs per-molecule identity?
- Does polymer handling preserve source-defined monomer or repeat-unit identity?
- Does any analysis or documentation rely on residue numbers alone where a
  chain/residue pair would be safer?
- Would a visualization user be able to map an interesting analysis result back
  to a specific unit in the generated topology?

The checklist is intentionally conceptual. The exact inspection method may vary
by topology format and analysis tool, but the invariants should remain the same.

## Contributor guidance

When contributing code or documentation that touches residue or chain identity,
preserve these invariants:

- Treat chain IDs as semantic component labels, not arbitrary formatting.
- Treat residues as chemically meaningful units for selection and aggregation.
- Prefer chain/residue pairs when referring to a specific unit, especially in
  systems with multiple chains or numbering overflow.
- Do not collapse distinct monomers, solvent molecules, ions, or co-solvents
  into a shared residue identity for convenience.
- Do not assume polymer residue numbering restarts per chain unless the source
  representation explicitly defines that behavior.
- Document public generated-topology behavior rather than private builder
  methods or internal implementation details.

These rules help new analyses remain compatible with the rest of PolyzyMD's
tooling and keep generated systems interpretable across simulation,
visualization, and comparison workflows.

## Summary

PolyzyMD's residue assignment convention preserves chemically meaningful
identity in generated topology files. The central idea is that one repeat unit is
one residue, interpreted as the smallest useful unit for selection, tracking, and
summarization.

The chain convention gives each major component a predictable role: chain A for
protein, chain B for substrate or ligand, chain C for polymer, and chains D and
later for solvent, ions, and remaining molecules. Together, chain IDs and
residue numbers provide stable handles for the units that analyses and
visualization tools need to reason about.
