# molplotly Cleanup & Modernization Plan

## Summary

molplotly is a single-module package (~520 lines) that adds interactive molecule hover tooltips to plotly figures using RDKit and Dash. The core idea is solid but the package is broken on PyPI, relies on a deprecated dependency (`JupyterDash`), has fragile packaging, minimal tests, and several unaddressed user issues.

---

## Phase 1: Critical Fixes (Get it working again)

### 1.1 Replace `JupyterDash` with `Dash`
- **Why**: `JupyterDash` is deprecated and causes doubled plots + port errors on Dash >= 2.10 (issue #31). This is the #1 breakage.
- **What**: Adopt the approach from PR #32 (open since Dec 2023, never reviewed) — replace `jupyter-dash` imports with standard `dash.Dash`. Modern Dash (>=2.11) natively supports Jupyter without `JupyterDash`.
- **Files**: `molplotly/main.py`, `tests/test_add_molecules.py`

### 1.2 Modernize packaging — migrate to `pyproject.toml`
- **Why**: `setup.py` is legacy, `setup_pip.py` uses the removed `distutils`, and the current PyPI release (1.1.8) has mismatched version metadata making it uninstallable (issue #35).
- **What**:
  - Create `pyproject.toml` with all metadata, dependencies, and build config
  - Delete `setup.py` and `setup_pip.py`
  - Set a correct, bumped version (e.g. 2.0.0 given the breaking `JupyterDash` removal)
  - Pin minimum dependency versions sensibly: `dash>=2.11.0`, `plotly>=5.0.0`, `rdkit`, `pandas`
  - Remove `jupyter-dash`, `werkzeug`, `ipykernel`, `nbformat` from required deps (no longer needed without JupyterDash; Dash pulls in werkzeug transitively)
- **Files**: new `pyproject.toml`, delete `setup.py`, delete `setup_pip.py`

### 1.3 Update CI/CD
- **Why**: CI uses Python 3.8 (EOL), `rdkit-pypi` (renamed), and `actions/setup-python@v2` (outdated).
- **What**:
  - Bump to Python 3.10+ (or test matrix 3.10/3.11/3.12)
  - Use `actions/checkout@v4`, `actions/setup-python@v5`
  - Remove miniconda setup (rdkit installs fine via pip now)
  - Install via `pip install .[test]` (which will use pyproject.toml)
- **Files**: `.github/workflows/test.yml`

---

## Phase 2: Code Quality

### 2.1 Clean up `main.py`
- **Bare except** (line ~369): Catch specific `Exception` or `ValueError` instead of bare `except:`
- **Unused imports / dead code**: Audit and remove
- **Type hints**: Add missing type hints to `test_groups`, `find_correct_column_order`, `find_grouping`
- **Docstrings**: Add/improve docstrings for the public `add_molecules` function and helpers
- **Files**: `molplotly/main.py`

### 2.2 Improve `__init__.py`
- **Why**: `from .main import *` exports everything including internal helpers
- **What**: Define `__all__` to only export `add_molecules` (the public API), or use explicit imports
- **Files**: `molplotly/__init__.py`

### 2.3 Add `__version__`
- Use `importlib.metadata` to expose `__version__` from the installed package metadata (single source of truth from `pyproject.toml`)
- **Files**: `molplotly/__init__.py`

---

## Phase 3: Testing

### 3.1 Expand test coverage
- **Why**: Currently 1 test that only checks `isinstance(app, JupyterDash)` — no functional validation
- **What**: Add tests for:
  - Basic scatter plot with molecules
  - Color column grouping
  - Symbol column grouping
  - Facet column support
  - Multiple SMILES columns
  - Reaction SMILES drawing
  - Edge cases: missing SMILES, invalid SMILES, empty dataframe
  - Caption formatting functions
  - The `find_grouping` logic (unit test the helper directly)
- **Files**: `tests/test_add_molecules.py` (or split into multiple test files)

---

## Phase 4: Documentation & Metadata

### 4.1 Update README
- Update installation instructions
- Note the `JupyterDash` → `Dash` migration (breaking change for users who import JupyterDash type)
- Update the "known issues" section (remove stale items, add current limitations)
- Update badges if any

### 4.2 Update CITATION.cff
- **Why**: Shows version 1.1.1 and date 2022-03-01, both stale
- **What**: Bump version and date to match the new release

### 4.3 Update example notebooks
- Verify notebooks still run with the updated code
- Remove/fix any `JupyterDash`-specific patterns

---

## Phase 5: Address Open Issues (nice-to-haves / future work)

These are lower priority but worth tracking:

| Issue | Description | Effort |
|-------|-------------|--------|
| #34 | Streamlit integration | Medium — would need a different rendering approach |
| #29 | Embed in existing Dash app | Medium — return a Dash component instead of a full app |
| #26 | Support for stacked bar charts / graph_objects | Small-Medium |
| #4 | Export as standalone HTML | Hard — fundamental limitation of needing a Dash server |
| #20 | Usage question (resolved) | Can be closed |
| #21 | make_subplots support | Medium |

---

## Suggested Execution Order

1. **Phase 1.1** — Replace JupyterDash (unblocks everything)
2. **Phase 1.2** — pyproject.toml migration (fixes packaging)
3. **Phase 2.1-2.3** — Code cleanup (while we're in the code)
4. **Phase 1.3** — Update CI (validates the above changes)
5. **Phase 3** — Expand tests
6. **Phase 4** — Docs & metadata
7. **Phase 5** — Feature work (separate PRs)

Phases 1-4 could reasonably be done as a single "v2.0.0" release given the breaking JupyterDash change.
