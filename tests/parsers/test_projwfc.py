# -*- coding: utf-8 -*-
# pylint: disable=unused-argument,redefined-outer-name
"""Tests for the `ProjwfcParser`."""
from __future__ import absolute_import

import os
import pytest

from aiida_quantumespresso.parsers.parse_raw.projwfc import parse_lowdin_charges


@pytest.mark.parametrize('test_file', ('default', 'spin'))
def test_parse_lowdin_charges(test_file, parser_fixture_path, data_regression):
    path = os.path.join(parser_fixture_path, 'projwfc', test_file, 'aiida.out')
    with open(path) as handle:
        data, spill_param = parse_lowdin_charges(handle.read().splitlines())
    data['spill'] = spill_param
    data_regression.check(data)
