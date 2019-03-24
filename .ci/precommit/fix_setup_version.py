"""Pre-commit script to ensure that version numbers in `setup.json` and `aiida_quantumespresso/__init__.py` match."""
# pylint: disable=wrong-import-position
import os
import json
import sys

FILEPATH = os.path.split(os.path.realpath(__file__))[0]

# Get current JSON content
FILEPATH_SETUP_JSON = os.path.join(FILEPATH, os.pardir, os.pardir, 'setup.json')
with open(FILEPATH_SETUP_JSON) as f:
    CONTENT_SETUP_JSON = json.load(f)

# Retrieve version from python package
sys.path.insert(0, os.path.join(FILEPATH, os.pardir, os.pardir))
import aiida_quantumespresso
VERSION = aiida_quantumespresso.__version__

CONTENT_SETUP_JSON['version'] = VERSION

# Rewrite JSON in a 'consistent' way (sorted, indented)
with open(FILEPATH_SETUP_JSON, 'w') as f:
    json.dump(CONTENT_SETUP_JSON, f, indent=4, sort_keys=True)
