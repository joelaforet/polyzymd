# Known Issues

## 1. Sphinx Incremental Build Limitations (Severity: Low, Documented)

**Symptom:** After adding a new page to a `toctree` directive, the sidebar
in built documentation doesn't show the new page (other pages appear stale).

**Root cause:** Sphinx incremental builds (`make html`) don't always detect
toctree structural changes. This is standard Sphinx behavior.

**Fix:** Always run `make clean html` instead of `make html` after adding
or removing toctree entries. This is documented in
`.opencode/instructions/documentation.md`.

---

## 2. GitHub Issue #20 — Analysis Module TODOs (Severity: Tracking)

**Symptom:** Various incomplete features and inconsistencies in the analysis
module.

**Details:** The refactor roadmap is kept by the maintainer outside the repository. Key items:
- Standardize analyzer inheritance
- Unify result formats
- Add comprehensive tests
- Fix bugs #1 and #2 above

---

## 3. Pre-existing LSP Type Errors (Severity: Low, Cosmetic)

**Symptom:** Pyright/Pylance reports many type errors in `config/schema.py`,
`builders/system_builder.py`, `simulation/runner.py`, and `cli/main.py`.

**Root cause:** These are mostly due to:
- Pydantic v2 `default_factory` type inference issues (false positives)
- OpenMM unit system lacking type stubs
- `Optional` vs runtime `None` handling patterns

**Impact:** These do NOT affect runtime behavior. The code runs correctly.
They are static analysis noise from missing type stubs for scientific packages.

**Fix approach:** Add `py.typed` marker and targeted `# type: ignore` comments,
or contribute type stubs for OpenMM/OpenFF. Low priority.

---

## Resolved Issues

### Config Hash Mismatch Warning (Resolved in v1.3.0)

**Was:** the plugin framework's cache validation printed "Config hash
mismatch" 66+ times instead of once.

**Resolution:** the plugin framework and its cache validation were removed.
The study API compares each stored result's `record.json`, config hash
included, once per replicate and measures the replicate again on a mismatch.

### Sphinx Doc Build Warnings (Resolved in v1.3.0)

**Was:** 196 Sphinx build warnings (195 duplicate object descriptions + 1
non-consecutive header level).

**Root cause:** `autodoc_typehints = "description"` combined with
`special-members: __init__` generated duplicate `attribute` directives for
Pydantic model fields and dataclass fields. Two docstrings also had reST
indentation issues.

**Resolution:** Added `autodoc-pydantic` extension for Pydantic v2 rendering.
Added `:no-index:` to `automodule` directives for all modules containing
Pydantic models or dataclasses. Fixed reST formatting in `PolymerBuilder` and
`EquilibrationStageConfig` docstrings. Fixed non-consecutive header level in
`cli_reference.md`. See `.opencode/instructions/documentation.md` for
prevention rules.
