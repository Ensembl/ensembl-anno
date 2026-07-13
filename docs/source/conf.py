# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import datetime
from importlib.metadata import version as pkg_version

sys.path.insert(0, os.path.abspath("../../src/python"))


# -- Project information -----------------------------------------------------
project = "ensembl-anno"
author = "ensembl@dev.org"
copyright_owner = "EMBL-European Bioinformatics Institute"
copyright_dates = "[2016-%d]" % datetime.datetime.now().year
copyright = copyright_dates + " " + copyright_owner
html_baseurl = "https://ensembl.github.io/ensembl-anno/"

try:
    release = pkg_version("ensembl-anno")
except Exception:
    release = "development"

version = ".".join(release.split(".")[:2])

# The master toctree document.
master_doc = "index"

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

html_logo = "_static/ensembl_anno_logo.png"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.duration",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
highlight_language = "bash"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

nitpicky = True

nitpick_ignore_regex = [
    ("py:class", r"pathlib\..*"),
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

# Defining autodoc functionality

autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "undoc-members": False,
    "show-inheritance": True,
    "special-members": "__init__, __call__",
}

autodoc_typehints = "description"
autodoc_typehints_format = "short"

autosectionlabel_prefix_document = True
autosummary_generate = True  # Turn on sphinx.ext.autosummary

# Set napoleon_use_param to True to format parameters as in the docstring
napoleon_use_param = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True
autosummary_imported_members = False


# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "navigation_depth": 3,
    "show_nav_level": 2,
    "show_toc_level": 2,
    "navbar_align": "left",
    "secondary_sidebar_items": [
        "page-toc",
        "edit-this-page",
    ],
}

html_theme_options = {
    "logo": {
        "image_light": "_static/ensembl_anno_logo.png",
        "image_dark": "_static/ensembl_anno_logo.png",
        "text": "Ensembl Anno",
    },
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_css_files = [
    "custom.css",
]


suppress_warnings = [
    "docstring",
    "ref.citation",  # Example: Suppress warnings about citations
    "image.nonlocal_uri",  # Example: Suppress warnings about non-local images
    # Add other warnings to suppress as needed
]


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "ensembl-anno.tex", "Ensembl-anno Base Library Documentation", [author], "manual"),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "ensembl-anno", "Ensembl-anno Base Library Documentation", [author], 1)]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "ensembl-anno",
        "Ensembl-anno Library Documentation",
        author,
        "ensembl-anno",
        "Ensembl-anno Base Library.",
        "Miscellaneous",
    ),
]
