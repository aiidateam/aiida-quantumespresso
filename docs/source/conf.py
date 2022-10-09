# -*- coding: utf-8 -*-
#
# aiida-wannier90 documentation build configuration file, created by
# sphinx-quickstart on Fri Oct 10 02:14:52 2014.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.
import time

# Load the dummy profile even if we are running locally, this way the documentation will succeed even if the current
# default profile of the AiiDA installation does not use a Django backend.
from aiida.manage.configuration import load_documentation_profile

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import aiida_quantumespresso

load_documentation_profile()

# -- Project information -----------------------------------------------------

project = 'aiida-quantumespresso'
copyright = f"""2014-{time.localtime().tm_year}, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of
Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)),
Switzerland. All rights reserved"""

# The full version, including alpha/beta/rc tags.
release = aiida_quantumespresso.__version__
# The short X.Y version.
version = '.'.join(aiida_quantumespresso.__version__.split('.')[:2])

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    'sphinx_click',
    'sphinx_design',
]

# Setting the intersphinx mapping to other readthedocs
intersphinx_mapping = {
    'python': ('https://docs.python.org/3.8', None),
    'aiida': ('http://aiida_core.readthedocs.io/en/latest/', None),
}

# Settings for the `sphinx_copybutton` extension
copybutton_selector = 'div:not(.no-copy)>div.highlight pre'
copybutton_prompt_text = r'>>> |\.\.\. |(?:\(.*\) )?\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: '
copybutton_prompt_is_regexp = True

nitpick_ignore = [('py:obj', 'module')]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    'github_url': 'https://github.com/aiidateam/aiida-quantumespresso',
    'twitter_url': 'https://twitter.com/aiidateam',
    'use_edit_page_button': True,
}
html_static_path = ['_static']
html_context = {
    'github_user': 'aiidateam',
    'github_repo': 'aiida-quantumespresso',
    'github_version': 'main',
    'doc_path': 'docs/source',
}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = 'images/logo_docs.png'
html_static_path = ['_static']
html_css_files = ['aiida-custom.css', 'aiida-qe-custom.css']

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
html_use_opensearch = 'http://aiida-quantumespresso.readthedocs.io'

# Language to be used for generating the HTML full-text search index.
# Sphinx supports the following languages:
#   'da', 'de', 'en', 'es', 'fi', 'fr', 'hu', 'it', 'ja'
#   'nl', 'no', 'pt', 'ro', 'ru', 'sv', 'tr'
html_search_language = 'en'

# Output file base name for HTML help builder.
htmlhelp_basename = 'aiida-quantumespressodoc'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #'preamble': '',

    # Latex figure (float) alignment
    #'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
# latex_documents = [
# ]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
# man_pages = [
# ]

# If true, show URL addresses after external links.
#man_show_urls = False

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
# texinfo_documents = [
# ]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#texinfo_no_detailmenu = False

# Warnings to ignore when using the -n (nitpicky) option
# We should ignore any python built-in exception, for instance
nitpick_ignore = [
    ('py:exc', 'ArithmeticError'),
    ('py:exc', 'AssertionError'),
    ('py:exc', 'AttributeError'),
    ('py:exc', 'BaseException'),
    ('py:exc', 'BufferError'),
    ('py:exc', 'DeprecationWarning'),
    ('py:exc', 'EOFError'),
    ('py:exc', 'EnvironmentError'),
    ('py:exc', 'Exception'),
    ('py:exc', 'FloatingPointError'),
    ('py:exc', 'FutureWarning'),
    ('py:exc', 'GeneratorExit'),
    ('py:exc', 'IOError'),
    ('py:exc', 'ImportError'),
    ('py:exc', 'ImportWarning'),
    ('py:exc', 'IndentationError'),
    ('py:exc', 'IndexError'),
    ('py:exc', 'KeyError'),
    ('py:exc', 'KeyboardInterrupt'),
    ('py:exc', 'LookupError'),
    ('py:exc', 'MemoryError'),
    ('py:exc', 'NameError'),
    ('py:exc', 'NotImplementedError'),
    ('py:exc', 'OSError'),
    ('py:exc', 'OverflowError'),
    ('py:exc', 'PendingDeprecationWarning'),
    ('py:exc', 'ReferenceError'),
    ('py:exc', 'RuntimeError'),
    ('py:exc', 'RuntimeWarning'),
    ('py:exc', 'StandardError'),
    ('py:exc', 'StopIteration'),
    ('py:exc', 'SyntaxError'),
    ('py:exc', 'SyntaxWarning'),
    ('py:exc', 'SystemError'),
    ('py:exc', 'SystemExit'),
    ('py:exc', 'TabError'),
    ('py:exc', 'TypeError'),
    ('py:exc', 'UnboundLocalError'),
    ('py:exc', 'UnicodeDecodeError'),
    ('py:exc', 'UnicodeEncodeError'),
    ('py:exc', 'UnicodeError'),
    ('py:exc', 'UnicodeTranslateError'),
    ('py:exc', 'UnicodeWarning'),
    ('py:exc', 'UserWarning'),
    ('py:exc', 'VMSError'),
    ('py:exc', 'ValueError'),
    ('py:exc', 'Warning'),
    ('py:exc', 'WindowsError'),
    ('py:exc', 'ZeroDivisionError'),
    ('py:obj', 'str'),
    ('py:obj', 'list'),
    ('py:obj', 'tuple'),
    ('py:obj', 'int'),
    ('py:obj', 'float'),
    ('py:obj', 'bool'),
    ('py:obj', 'Mapping'),
]
