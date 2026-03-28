# Known Issues

## 1. Config Hash Mismatch Warning (Severity: Medium)

**Symptom:** When running analysis, the warning "Config hash mismatch" prints
66+ times to the console instead of once.

**Root cause:** The cache validation logic checks hashes per-frame or per-chunk
rather than once at the start. Each check triggers the warning independently.

**Location:** `src/polyzymd/analyses/shared/config_hash.py` (cache validation logic)

**Fix approach:** Add a "warned once" flag or move the hash check to a
single-call validation step before the analysis loop begins.

**Workaround:** Ignore the warning output; the analysis results are still
correct.

---

## 2. Contact Criteria Cutoff Mismatch (Severity: Medium)

**Symptom:** Cached analysis results may use a 4.0A cutoff while the current
config specifies 4.5A (or vice versa). This leads to incorrect results being
loaded from cache.

**Root cause:** The cache key does not include the cutoff value, so changing
the cutoff doesn't invalidate the cache.

**Location:** `src/polyzymd/analyses/contacts/` and `analyses/shared/config_hash.py`

**Fix approach:** Include the cutoff value (and other analysis parameters)
in the cache key hash computation.

**Workaround:** Manually clear the cache directory before re-running with
different cutoff values.

---

## 3. Docs Sidebar Stale After Adding Pages (Severity: Low)

**Symptom:** After adding a new page to a `toctree` directive, the sidebar
in built documentation doesn't show the new page (other pages appear stale).

**Root cause:** Sphinx incremental builds (`make html`) don't always detect
toctree structural changes.

**Fix:** Always run `make clean html` instead of `make html` after adding
or removing toctree entries.

---

## 4. GitHub Issue #20 — Analysis Module TODOs (Severity: Tracking)

**Symptom:** Various incomplete features and inconsistencies in the analysis
module.

**Details:** See `analysis-module.md` for the full roadmap. Key items:
- Standardize analyzer inheritance
- Unify result formats
- Add comprehensive tests
- Fix bugs #1 and #2 above

---

## 5. Pre-existing LSP Type Errors (Severity: Low, Cosmetic)

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
