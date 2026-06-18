# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Cardillo Documentation"
copyright = "2026, Jan Kautz, Jonas Breuling"
author = "Jan Kautz, Jonas Breuling"
release = "1.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import subprocess

# Auto-run Doxygen before every Sphinx build
subprocess.call("cd .. && doxygen docs/Doxyfile", shell=True)

extensions = [
    "breathe",
    "exhale",
    "sphinx.ext.mathjax",  # LaTeX math support
    "sphinx.ext.graphviz",  # call graphs
    "sphinxcontrib.mermaid",  # Mermaid flow diagrams
]

html_theme = "sphinx_rtd_theme"

# --- Breathe ---
breathe_projects = {"Cardillo": "./_doxygen/xml"}
breathe_default_project = "Cardillo"

# --- Exhale ---
exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "API Reference",
    "doxygenStripFromPath": "..",
    "createTreeView": True,
    "listingExclude": [
        "function",
        "struct",
        "variable",
        "define",
        "enum",
        "typedef",
        "file",
        "dir",
        "union",
    ],
    "contentsDirectives": False,
    "generateBreatheFileDirectives": False,
}

# --- Math ---
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"

# --- Mermaid ---
mermaid_version = "11.4.1"
