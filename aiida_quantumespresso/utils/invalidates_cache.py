#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compatibility module for using the 'invalidates_cache' functionality.

Helper to allow using the 'invalidates_cache' functionality without
breaking compatibility with AiiDA 1.0. To mark an exit code as
invalidating cache, pass **INVALIDATES_CACHE to the spec.exit_code call.
If the functionality is not implemented in the used AiiDA version, this
has no effect.
This module should be removed once AiiDA 1.0 (and python2) compatibility
is no longer maintained, instead passing 'invalidates_cache=True'
directly.
"""

from __future__ import absolute_import
from aiida.engine.processes.exit_code import ExitCode

__all__ = ('INVALIDATES_CACHE_INPUT',)

if 'invalidates_cache' in ExitCode._fields:
    INVALIDATES_CACHE_INPUT = {'invalidates_cache': True}
else:
    INVALIDATES_CACHE_INPUT = {}
