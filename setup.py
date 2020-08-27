# -*- coding: utf-8 -*-
"""Setup for the `aiida-quantumespresso` plugin which provides an interface to Quantum ESPRESSO for `aiida-core`."""
try:
    import fastentrypoints  # pylint: disable=unused-import
except ImportError:
    # This should only occur when building the package, i.e. for `python setup.py sdist/bdist_wheel`
    pass


def setup_package():
    """Set up the package."""
    import json
    from setuptools import setup, find_packages

    filename_setup_json = 'setup.json'
    filename_description = 'README.md'

    with open(filename_setup_json, 'r') as handle:
        setup_json = json.load(handle)

    with open(filename_description, 'r') as handle:
        description = handle.read()

    setup(
        include_package_data=True,
        packages=find_packages(),
        reentry_register=True,
        long_description=description,
        long_description_content_type='text/markdown',
        **setup_json
    )


if __name__ == '__main__':
    setup_package()
