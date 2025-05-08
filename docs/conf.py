import os
import sys

sys.path.insert(0, os.path.abspath('.'))

try:
    import picht
    release = picht.__version__
except ImportError:
    release = '1.1.4'

project = 'picht'
copyright = '2025, Rishiit Sharma'
author = 'Rishiit Sharma'
version = '.'.join(release.split('.')[:2])

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
]

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False,
    'special-members': False,
    'inherited-members': False,
    'show-inheritance': True,
}

autosummary_generate = True