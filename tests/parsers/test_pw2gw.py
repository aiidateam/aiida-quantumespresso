# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,redefined-outer-name
"""Tests for the `Pw2gwParser`."""
from __future__ import absolute_import

from aiida import orm


def test_pw2gw_default(
    aiida_profile,
    fixture_localhost,
    generate_parser,
    generate_calc_job_node,
):
    """Test a normal pw2gw.x output."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.pw2gw'
    entry_point_parser = 'quantumespresso.pw2gw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name)

    parser = generate_parser(entry_point_parser)

    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)


def test_pw2gw_failed_missing(
    aiida_profile,
    fixture_localhost,
    generate_parser,
    generate_calc_job_node,
):
    """Test a pw2gw.x output where file are missing."""
    name = 'failed_missing'
    entry_point_calc_job = 'quantumespresso.pw2gw'
    entry_point_parser = 'quantumespresso.pw2gw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name)

    parser = generate_parser(entry_point_parser)

    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status, node.process_class.exit_codes.ERROR_OUTPUT_FILES.status
    assert orm.Log.objects.get_logs_for(node)


def test_pw2gw_failed_corrupted_file(
    aiida_profile,
    fixture_localhost,
    generate_parser,
    generate_calc_job_node,
):
    """Test a pw2gw.x output where file are corrupted."""
    name = 'failed_corrupted_file'
    entry_point_calc_job = 'quantumespresso.pw2gw'
    entry_point_parser = 'quantumespresso.pw2gw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name)

    parser = generate_parser(entry_point_parser)

    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status, node.process_class.exit_codes.ERROR_OUTPUT_FILES.status
    assert orm.Log.objects.get_logs_for(node)
