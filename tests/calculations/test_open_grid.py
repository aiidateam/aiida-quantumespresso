"""Tests for the `OpenGridCalculation` class."""

import pytest
from aiida import orm
from aiida.common import AttributeDict, datastructures
from aiida.common.exceptions import InputValidationError


@pytest.fixture
def generate_inputs(fixture_localhost, fixture_sandbox, fixture_code, generate_remote_data):
    """Return only those inputs that the parser will expect to be there."""

    def _generate_inputs(parameters=None, settings=None):
        from aiida_quantumespresso.utils.resources import get_default_options

        if parameters is None:
            parameters = {'INPUTPP': {}}

        return AttributeDict(
            {
                'code': fixture_code('quantumespresso.open_grid'),
                'parent_folder': generate_remote_data(
                    fixture_localhost, fixture_sandbox.abspath, 'quantumespresso.open_grid'
                ),
                'parameters': orm.Dict(dict=parameters),
                'settings': orm.Dict(dict=settings),
                'metadata': {'options': get_default_options()},
            }
        )

    return _generate_inputs


def test_open_grid_default(fixture_sandbox, generate_calc_job, generate_inputs, file_regression):
    """Test a default `OpenGridCalculation`."""
    entry_point_name = 'quantumespresso.open_grid'
    inputs = generate_inputs()

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['aiida.out']

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


@pytest.mark.parametrize(
    ('parameters', 'message'),
    [
        (
            {'INPUTPP': {'outdir': './out/'}},
            r"You cannot specify explicitly the 'outdir' key in the 'INPUTPP' namelist",
        ),
        (
            {'INPUTPP': {'overwrite_prefix': True}},
            r"You cannot specify explicitly the 'overwrite_prefix' key in the 'INPUTPP' namelist",
        ),
        ({'INPUTPP': {'prefix': 'aiida'}}, r"You cannot specify explicitly the 'prefix' key in the 'INPUTPP' namelist"),
    ],
)
def test_open_grid_invalid_parameters(fixture_sandbox, generate_calc_job, generate_inputs, parameters, message):
    """Test that launching `OpenGridCalculation` fails for invalid parameters."""
    entry_point_name = 'quantumespresso.open_grid'

    inputs = generate_inputs(parameters=parameters)

    with pytest.raises(InputValidationError, match=message) as exception:
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    assert message in str(exception.value)
