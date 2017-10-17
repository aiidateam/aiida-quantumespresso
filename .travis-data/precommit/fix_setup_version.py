import os
import json
import sys

this_path = os.path.split(os.path.realpath(__file__))[0]

# Get current JSON content
setup_path = os.path.join(this_path, os.pardir, os.pardir, 'setup.json')
with open(setup_path) as f:
	setup_content = json.load(f)

# Retrieve version from python package
sys.path.insert(0, os.path.join(this_path, os.pardir, os.pardir))
import aiida_quantumespresso
version = aiida_quantumespresso.__version__

setup_content['version'] = version

# Rewrite JSON in a 'consistent' way (sorted, indented)
with open(setup_path, 'w') as f:
	json.dump(setup_content, f, indent=2, sort_keys=True)
