# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `PpCalculation` class."""
import pytest

from aiida import orm
from aiida.common import datastructures, AttributeDict


@pytest.fixture
def generate_inputs(fixture_localhost, fixture_sandbox, fixture_code, generate_remote_data):
    """Return only those inputs that the parser will expect to be there."""

    def _generate_inputs(parameters=None, settings=None):
        from aiida_quantumespresso.utils.resources import get_default_options

        if parameters is None:
            parameters = {'INPUTPP': {'plot_num': 1}, 'PLOT': {'iflag': 3}}

        return AttributeDict({
            'code':
            fixture_code('quantumespresso.pp'),
            'parent_folder':
            generate_remote_data(fixture_localhost, fixture_sandbox.abspath, 'quantumespresso.pw'),
            'parameters':
            orm.Dict(dict=parameters),
            'settings':
            orm.Dict(dict=settings),
            'metadata': {
                'options': get_default_options()
            }
        })

    return _generate_inputs


def test_pp_default(fixture_sandbox, generate_calc_job, generate_inputs, file_regression):
    """Test a default `PpCalculation`."""
    entry_point_name = 'quantumespresso.pp'
    inputs = generate_inputs()

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['aiida.out']
    retrieve_temporary_list = ['aiida.fileout', ('aiida.filplot_*aiida.fileout', '.', 0)]
    local_copy_list = []

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert len(calc_info.retrieve_temporary_list) == 2
    for element in retrieve_temporary_list:
        assert element in calc_info.retrieve_temporary_list

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_pp_keep_plot_file(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test a `PpCalculation` where we want to retrieve the plot file."""
    entry_point_name = 'quantumespresso.pp'
    inputs = generate_inputs()
    inputs.metadata.options.keep_plot_file = True

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    retrieve_list = ['aiida.out', 'aiida.fileout', ('aiida.filplot_*aiida.fileout', '.', 0)]
    retrieve_temporary_list = []
    local_copy_list = []

    # Check the attributes of the returned `CalcInfo`, no need to check the input file as it is not affected
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)
    assert len(calc_info.retrieve_list) == 3
    for element in retrieve_list:
        assert element in calc_info.retrieve_list


def test_pp_cmdline_setting(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test a `PpCalculation` with user-defined cmdline settings."""
    entry_point_name = 'quantumespresso.pp'
    inputs = generate_inputs(settings={'cmdline': ['-npools', '2']})
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    assert calc_info.codes_info[0].cmdline_params == ['-npools', '2']


@pytest.mark.parametrize(
    ('parameters', 'message'),
    (
        ({}, 'parameter `INPUTPP.plot_num` must be explicitly set'),
        ({
            'INPUTPP': {}
        }, 'parameter `INPUTPP.plot_num` must be explicitly set'),
        ({
            'INPUTPP': {
                'plot_num': 'str'
            }
        }, '`INTPUTPP.plot_num` must be an integer in the range'),
        ({
            'INPUTPP': {
                'plot_num': 14
            }
        }, '`INTPUTPP.plot_num` must be an integer in the range'),
        ({
            'INPUTPP': {
                'plot_num': 1
            }
        }, 'parameter `PLOT.iflag` must be explicitly set'),
        ({
            'INPUTPP': {
                'plot_num': 1
            },
            'PLOT': {}
        }, 'parameter `PLOT.iflag` must be explicitly set'),
        ({
            'INPUTPP': {
                'plot_num': 1
            },
            'PLOT': {
                'iflag': 'str'
            }
        }, '`PLOT.iflag` must be an integer in the range 0-4'),
        ({
            'INPUTPP': {
                'plot_num': 1
            },
            'PLOT': {
                'iflag': 5
            }
        }, '`PLOT.iflag` must be an integer in the range 0-4'),
    ),
)
def test_pp_invalid_parameters(fixture_sandbox, generate_calc_job, generate_inputs, parameters, message):
    """Test that launching `PpCalculation` fails for invalid parameters."""
    entry_point_name = 'quantumespresso.pp'

    with pytest.raises(ValueError) as exception:
        inputs = generate_inputs(parameters=parameters)
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    assert message in str(exception.value)
