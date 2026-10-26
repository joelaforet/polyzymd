# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the source directory to the path so autodoc can find modules
sys.path.insert(0, os.path.abspath("../../src"))

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
    "mbuild",
]

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "PolyzyMD"
copyright = "2026, Joseph R. Laforet Jr."
author = "Joseph R. Laforet Jr."
release = "1.3.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # Auto-generate docs from docstrings
    "sphinxcontrib.autodoc_pydantic",
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
napoleon_use_ivar = False
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

# autodoc-pydantic settings (Pydantic v2 model rendering)
autodoc_pydantic_model_show_field_summary = False
autodoc_pydantic_model_show_json = False
autodoc_pydantic_model_show_config_summary = False
autodoc_pydantic_model_show_config_member = False
autodoc_pydantic_model_show_validator_summary = True
autodoc_pydantic_model_show_validator_members = True
autodoc_pydantic_field_show_constraints = True
autodoc_pydantic_field_show_required = True
autodoc_pydantic_field_show_default = True
autodoc_pydantic_field_show_alias = False
autodoc_pydantic_model_undoc_members = True

# Autosummary settings
autosummary_generate = True

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

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_extra_path = ["_extra"]
html_title = "PolyzyMD"
html_favicon = "_static/favicon.png"
# The badge has a transparent background, so it sits on the header in light and dark mode.
html_logo = "_static/polyzymd_badge_256.png"
html_theme_options = {
    "logo": {"text": "PolyzyMD", "alt_text": "PolyzyMD home"},
    # The top-level toctree (Get started, Tutorials, How-to, Explanation,
    # Reference, Contributing) becomes the header tabs: one per Diataxis section.
    "header_links_before_dropdown": 6,
    "navbar_align": "left",
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/joelaforet/polyzymd",
            "icon": "fa-brands fa-github",
        },
    ],
    "show_toc_level": 2,
    "navigation_with_keys": False,
    "show_prev_next": True,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "pygments_light_style": "tango",
    "pygments_dark_style": "monokai",
}
# The landing page has no left sidebar: its cards are the navigation.
html_sidebars = {"index": []}
html_show_sourcelink = False
html_css_files = [
    "custom.css",
]
# Looping molecular movies pause for readers who ask for reduced motion.
html_js_files = [
    "pz-movie.js",
]
# Rank task pages (get started, tutorials, how-to) above API pages, and pages
# with the query words in their title above pages that mention them once.
# scripts/check_docs_search.py measures the ranking on the queries users type.
html_search_scorer = os.path.join(os.path.dirname(os.path.abspath(__file__)), "search_scorer.js")

# -- Options for TODO extension ----------------------------------------------
todo_include_todos = True

# -- Options for copy button -------------------------------------------------
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
