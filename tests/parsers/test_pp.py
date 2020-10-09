# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,redefined-outer-name
"""Tests for the `PpParser`."""
import pytest

from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def generate_inputs_1d():
    """Return inputs that parser will expect for the 1D case."""

    inputs = {
        'parent_folder':
        orm.FolderData().store(),
        'parameters':
        orm.Dict(
            dict={
                'INPUTPP': {
                    'plot_num': 11
                },
                'PLOT': {
                    'iflag': 1,
                    'e1': [[1, 1.0], [2, 0.0], [3, 0.0]],
                    'x0': [[1, 0.], [2, 0.], [3, 0.]],
                    'nx': 100
                }
            }
        )
    }
    return AttributeDict(inputs)


@pytest.fixture
def generate_inputs_1d_spherical():
    """Return inputs that parser will expect for the 1D case with spherical averaging."""

    inputs = {
        'parent_folder':
        orm.FolderData().store(),
        'parameters':
        orm.Dict(
            dict={
                'INPUTPP': {
                    'plot_num': 11
                },
                'PLOT': {
                    'iflag': 0,
                    'e1': [[1, 1.0], [2, 0.0], [3, 0.0]],
                    'x0': [[1, 0.], [2, 0.], [3, 0.]],
                    'nx': 100
                }
            }
        )
    }
    return AttributeDict(inputs)


@pytest.fixture
def generate_inputs_2d():
    """Return inputs that parser will expect for the 2D case."""

    inputs = {
        'parent_folder':
        orm.FolderData().store(),
        'parameters':
        orm.Dict(
            dict={
                'INPUTPP': {
                    'plot_num': 11
                },
                'PLOT': {
                    'iflag': 2,
                    'e1': [[1, 1.0], [2, 0.0], [3, 0.0]],
                    'e2': [[1, 0.0], [2, 1.0], [3, 0.0]],
                    'x0': [[1, 0.], [2, 0.], [3, 0.]],
                    'nx': 10,
                    'ny': 10
                }
            }
        )
    }
    return AttributeDict(inputs)


@pytest.fixture
def generate_inputs_polar():
    """Return inputs that parser will expect for the polar case."""

    inputs = {
        'parent_folder':
        orm.FolderData().store(),
        'parameters':
        orm.Dict(
            dict={
                'INPUTPP': {
                    'plot_num': 11
                },
                'PLOT': {
                    'iflag': 4,
                    'e1': [[1, 1.0], [2, 0.0], [3, 0.0]],
                    'x0': [[1, 0.], [2, 0.], [3, 0.]],
                    'nx': 100
                }
            }
        )
    }
    return AttributeDict(inputs)


@pytest.fixture
def generate_inputs_3d():
    """Return inputs that parser will expect for the 3D case."""

    inputs = {
        'parent_folder': orm.FolderData().store(),
        'parameters': orm.Dict(dict={
            'INPUTPP': {
                'plot_num': 11
            },
            'PLOT': {
                'iflag': 3,
            }
        })
    }
    return AttributeDict(inputs)


def test_pp_default_1d(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_1d, data_regression, num_regression
):
    """Test a default `pp.x` calculation producing a 1D data set."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default_1d', generate_inputs_1d)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_data' in results
    assert len(results['output_data'].get_arraynames()) == 4
    data_array = results['output_data'].get_array('data')
    coords_array = results['output_data'].get_array('x_coordinates')
    units_array = results['output_data'].get_array('x_coordinates_units')
    num_regression.check({
        'data_array': data_array,
        'coords_array': coords_array
    },
                         default_tolerance=dict(atol=0, rtol=1e-18))
    data_regression.check({'parameters': results['output_parameters'].get_dict(), 'units_array': units_array.tolist()})


def test_pp_default_1d_spherical(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_1d_spherical, data_regression,
    num_regression
):
    """Test a default `pp.x` calculation producing a 1D data set with spherical averaging."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(
        entry_point_calc_job, fixture_localhost, 'default_1d_spherical', generate_inputs_1d_spherical
    )
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_data' in results
    assert len(results['output_data'].get_arraynames()) == 6
    data_array = results['output_data'].get_array('data')
    data_array_int = results['output_data'].get_array('integrated_data')
    coords_array = results['output_data'].get_array('x_coordinates')
    data_units_array = results['output_data'].get_array('data_units')
    data_int_units_array = results['output_data'].get_array('integrated_data_units')
    coords_units_array = results['output_data'].get_array('x_coordinates_units')
    num_regression.check({
        'data_array': data_array,
        'data_array_int': data_array_int,
        'coords_array': coords_array
    },
                         default_tolerance=dict(atol=0, rtol=1e-18))
    data_regression.check({
        'parameters': results['output_parameters'].get_dict(),
        'data_units_array': data_units_array.tolist(),
        'data_int_units_array': data_int_units_array.tolist(),
        'coords_units_array': coords_units_array.tolist(),
    })


def test_pp_default_2d(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_2d, data_regression, num_regression
):
    """Test a default `pp.x` calculation producing a 2D data set."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default_2d', generate_inputs_2d)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_data' in results
    assert len(results['output_data'].get_arraynames()) == 4
    data_array = results['output_data'].get_array('data').flatten()
    coords_array = results['output_data'].get_array('xy_coordinates').flatten()
    data_units_array = results['output_data'].get_array('data_units')
    coords_units_array = results['output_data'].get_array('xy_coordinates_units')
    num_regression.check({
        'data_array': data_array,
        'coords_array': coords_array
    },
                         default_tolerance=dict(atol=0, rtol=1e-18))
    data_regression.check({
        'parameters': results['output_parameters'].get_dict(),
        'data_units': data_units_array.tolist(),
        'coords_units': coords_units_array.tolist()
    })


def test_pp_default_polar(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_polar, data_regression, num_regression
):
    """Test a default `pp.x` calculation producing a polar coordinates data set."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default_polar', generate_inputs_polar)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_data' in results
    assert len(results['output_data'].get_arraynames()) == 2
    data_array = results['output_data'].get_array('data')
    data_units_array = results['output_data'].get_array('data_units')
    num_regression.check({
        'data_array': data_array,
    }, default_tolerance=dict(atol=0, rtol=1e-18))
    data_regression.check({
        'parameters': results['output_parameters'].get_dict(),
        'data_units': data_units_array.tolist(),
    })


def test_pp_default_3d(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_3d, data_regression, num_regression
):
    """Test a default `pp.x` calculation producing a 3D data set."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default_3d', generate_inputs_3d)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_data' in results
    assert len(results['output_data'].get_arraynames()) == 4
    data_array = results['output_data'].get_array('data').flatten()
    voxel_array = results['output_data'].get_array('voxel').flatten()
    data_units_array = results['output_data'].get_array('data_units')
    coordinates_units_array = results['output_data'].get_array('coordinates_units')
    num_regression.check({
        'data_array': data_array,
        'voxel_array': voxel_array,
    },
                         default_tolerance=dict(atol=0, rtol=1e-18))
    data_regression.check({
        'parameters': results['output_parameters'].get_dict(),
        'data_units': data_units_array.tolist(),
        'coordinates_units': coordinates_units_array.tolist()
    })


def test_pp_default_3d_keep_plot_file(generate_calc_job_node, generate_parser, generate_inputs_3d, tmpdir):
    """Test a `pp.x` calculation where `keep_plot_file=False` meaning files will be parsed from temporary directory."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    attributes = {'options': {'keep_plot_file': False}, 'retrieve_temporary_list': ['aiida.fileout']}
    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name='default_3d',
        inputs=generate_inputs_3d,
        attributes=attributes,
        retrieve_temporary=(tmpdir, ['aiida.fileout'])
    )
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_data' in results
    assert len(results['output_data'].get_arraynames()) == 4


def test_pp_default_3d_multiple(generate_calc_job_node, generate_parser, generate_inputs_3d):
    """Test a default `pp.x` calculation producing multiple files in 3D format."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(entry_point_calc_job, test_name='default_3d_multiple', inputs=generate_inputs_3d)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_data_multiple' in results

    # Since the actual parsing of the file content itself is done in `test_pp_default_3d` we will not do it here again
    for key in ['K001_B001', 'K001_B002']:
        assert key in results['output_data_multiple']
        node = results['output_data_multiple'][key]
        assert len(node.get_arraynames()) == 4


def test_pp_default_3d_failed_missing(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_3d):
    """Test a default `pp.x` calculation where no files are retrieved, or StdOut is missing."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(
        entry_point_calc_job, fixture_localhost, 'default_3d_failed_missing', generate_inputs_3d
    )
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    # This exception always fires first, although specifically it is for a missing StdOut file
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_MISSING.status


def test_pp_default_3d_failed_missing_data(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_3d
):
    """Test a default `pp.x` calculation where the aiida.fileout file is missing."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(
        entry_point_calc_job, fixture_localhost, 'default_3d_failed_missing_data', generate_inputs_3d
    )
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_DATAFILE_MISSING.status


def test_pp_default_3d_failed_interrupted(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_3d
):
    """Test a default `pp.x` calculation where the StdOut file is present but incomplete."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(
        entry_point_calc_job, fixture_localhost, 'default_3d_failed_interrupted', generate_inputs_3d
    )
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE.status


def test_pp_default_3d_failed_format(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs_3d):
    """Test a default `pp.x` calculation where an unsupported output file format is used."""
    entry_point_calc_job = 'quantumespresso.pp'
    entry_point_parser = 'quantumespresso.pp'

    node = generate_calc_job_node(
        entry_point_calc_job, fixture_localhost, 'default_3d_failed_format', generate_inputs_3d
    )
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_DATAFILE_PARSE.status
