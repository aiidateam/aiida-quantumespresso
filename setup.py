# -*- coding: utf-8 -*-
"""Setup script for aiida-core package."""
from __future__ import absolute_import
import json

# pylint: disable=wrong-import-order,unused-import
import fastentrypoints
from setuptools import setup, find_packages

if __name__ == '__main__':
    with open('setup.json', 'r') as info:
        SETUP_JSON = json.load(info)
    setup(include_package_data=True, reentry_register=True, packages=find_packages(), **SETUP_JSON)
