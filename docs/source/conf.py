# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import re
import sys

# Add the source directory to the path so autodoc can find modules
sys.path.insert(0, os.path.abspath("../../src"))


def _get_version() -> str:
    """Read the project version, preferring installed metadata with __init__.py fallback."""
    try:
        from importlib.metadata import version

        return version("polyzymd")
    except Exception:
        # Fallback: parse __init__.py directly (works without pip install)
        init = os.path.join(os.path.dirname(__file__), "..", "..", "src", "polyzymd", "__init__.py")
        with open(init) as f:
            match = re.search(r'__version__\s*=\s*"([^"]+)"', f.read())
        return match.group(1) if match else "unknown"


# Mock imports for packages that may not be available during doc build
# This allows autodoc to generate documentation without actually importing these
# heavy dependencies (OpenMM, OpenFF, etc.) which require CUDA/GPU libraries
autodoc_mock_imports = [
    "MDAnalysis",
    "MDAnalysis.analysis",
    "MDAnalysis.analysis.base",
    "MDAnalysis.lib",
    "openmm",
    "openff",
    "openff.toolkit",
    "openff.units",
    "openff.interchange",
    "openmmtools",
    "openmmforcefields",
    "pdbfixer",
    "mdtraj",
    "parmed",
    "rdkit",
    "numpy",
    "polymerist",
]

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "PolyzyMD"
copyright = "2026, Joseph R. Laforet Jr."
author = "Joseph R. Laforet Jr."
release = _get_version()
version = ".".join(release.split(".")[:2])  # short X.Y version for Sphinx

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # Auto-generate docs from docstrings
    "sphinx.ext.autosummary",  # Generate summary tables
    "sphinx.ext.napoleon",  # Support Google/NumPy style docstrings
    "sphinx.ext.viewcode",  # Add links to source code
    "sphinx.ext.intersphinx",  # Link to other projects' documentation
    "sphinx.ext.todo",  # Support TODO directives
    "sphinx.ext.mathjax",  # Render LaTeX math with MathJax
    "sphinx_copybutton",  # Add copy button to code blocks
    "myst_parser",  # Support Markdown files
    "sphinx_design",  # Tabs, cards, grids for better UX
]

# Napoleon settings for docstring parsing
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
    "show-inheritance": True,
}
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"

# Autosummary settings
# Stub generation is disabled because all API docs use automodule directives
# directly.  Enabling it would generate stub pages that duplicate every
# class/attribute already documented in the api/*.md pages.
autosummary_generate = False

# Intersphinx mapping to external documentation
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "openmm": ("http://docs.openmm.org/latest/api-python/", None),
}

# MyST parser settings (for Markdown support)
myst_enable_extensions = [
    "colon_fence",  # ::: fence for directives
    "deflist",  # Definition lists
    "fieldlist",  # Field lists
    "tasklist",  # Task lists with checkboxes
    "dollarmath",  # Enable $...$ and $$...$$ math syntax
]
myst_heading_anchors = 3

# Templates path
templates_path = ["_templates"]

# Patterns to exclude
exclude_patterns = []

# Source file suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Logo configuration
html_logo = "_static/logo.png"
html_favicon = "_static/logo.png"

# Theme options
html_theme_options = {
    "logo_only": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": True,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 3,
    "includehidden": True,
    "titles_only": False,
}

# Custom CSS
html_css_files = [
    "custom.css",
]

# -- Options for TODO extension ----------------------------------------------
todo_include_todos = True

# -- Options for copy button -------------------------------------------------
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
