import os
import sys
import matplotlib
matplotlib.use('agg')

try:
   import picht
   release = picht.__version__
except ImportError:
   release = '2.1.0'

version = '.'.join(release.split('.')[:2])

project = 'picht'
copyright = '2025, Rishiit Sharma'
author = 'Rishiit Sharma'

extensions = [
   'sphinx.ext.viewcode',
   'sphinx.ext.napoleon',
   'sphinx.ext.autosummary',
   'sphinx_gallery.gen_gallery',
]

napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_attr_annotations = True

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
   'special-members': '__init__',
   'inherited-members': False,
   'show-inheritance': True,
}

autosummary_generate = True

def custom_gallery_order(example_file):
   import os
   import re
   filename = os.path.basename(example_file)
   match = re.match(r'(\d+)_.*\.py$', filename)
   if match:
       return int(match.group(1))
   return float('inf')

sphinx_gallery_conf = {
   'examples_dirs': '../examples',
   'gallery_dirs': 'auto_examples',
   'filename_pattern': r'.*\.py$',
   'plot_gallery': True,
   'download_all_examples': True,
   'line_numbers': True,
   'remove_config_comments': True,
   'thumbnail_size': (400, 300),
   'min_reported_time': 0,
   'within_subsection_order': custom_gallery_order,
}

intersphinx_mapping = {
   'python': ('https://docs.python.org/3', None),
   'numpy': ('https://numpy.org/doc/stable/', None),
   'scipy': ('https://docs.scipy.org/doc/scipy/', None),
   'matplotlib': ('https://matplotlib.org/stable/', None),
}

import logging
logging.getLogger('sphinx_gallery').setLevel(logging.DEBUG)