# -*- coding: utf-8 -*-
"""Tests for the `MatdynParser`."""
from aiida import orm
from aiida.common import AttributeDict


def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    kpoints = orm.KpointsData()
    kpoints.set_kpoints_mesh([1, 1, 1])

    return AttributeDict({
        'kpoints': kpoints,
    })


def test_matdyn_default(fixture_localhost, generate_calc_job_node, generate_parser, data_regression):
    """Test a default `matdyn.x` calculation."""
    entry_point_calc_job = 'quantumespresso.matdyn'
    entry_point_parser = 'quantumespresso.matdyn'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_phonon_bands' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_phonon_bands': results['output_phonon_bands'].attributes
    })
