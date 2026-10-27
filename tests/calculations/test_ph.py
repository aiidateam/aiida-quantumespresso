"""Tests for the `PhCalculation` class."""

from pathlib import Path

import pytest
from aiida import orm
from aiida.common import datastructures
from aiida.plugins import CalculationFactory

PwCalculation = CalculationFactory('quantumespresso.pw')
PhCalculation = CalculationFactory('quantumespresso.ph')


def test_ph_default(fixture_sandbox, generate_inputs_ph, generate_calc_job, file_regression):
    """Test a default `PhCalculation`."""
    entry_point_name = 'quantumespresso.ph'
    inputs = generate_inputs_ph()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    cmdline_params = ['-in', 'aiida.in']
    retrieve_list = ['./out/_ph0/aiida.phsave/tensors.xml', 'DYN_MAT', 'aiida.out']
    local_copy_list = []

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.codes_info[0].cmdline_params) == sorted(cmdline_params)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.remote_symlink_list) == sorted([])

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['DYN_MAT', 'aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_ph_qpoint_list(
    fixture_sandbox,
    generate_inputs_ph,
    generate_calc_job,
    generate_structure,
    generate_kpoints_mesh,
    file_regression,
):
    """Test a `PhCalculation` with a qpoint list instead of a mesh."""
    entry_point_name = 'quantumespresso.ph'

    structure = generate_structure()
    kpoints = generate_kpoints_mesh(2).get_kpoints_mesh(print_list=True)
    qpoints = orm.KpointsData()
    qpoints.set_cell(structure.cell)
    qpoints.set_kpoints(kpoints)

    inputs = generate_inputs_ph()
    inputs['qpoints'] = qpoints
    generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_ph_initialization_only(fixture_sandbox, generate_inputs_ph, generate_calc_job):
    """Test a ``PhCalculation`` that should just run the initialization."""
    entry_point_name = 'quantumespresso.ph'
    inputs = generate_inputs_ph()
    inputs['settings'] = orm.Dict({'only_initialization': True})
    generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    assert (Path(fixture_sandbox.abspath) / f'{PhCalculation._PREFIX}.EXIT').exists()


@pytest.mark.parametrize('symlink', (True, False))
@pytest.mark.parametrize('electron_phonon', (None, 'interpolated'))
def test_ph_restart(
    fixture_sandbox,
    fixture_localhost,
    generate_inputs_ph,
    generate_calc_job,
    generate_remote_data,
    tmp_path,
    symlink,
    electron_phonon,
):
    """Test a ``PhCalculation`` that restarts from the ``parent_folder`` of another ``PhCalculation``.

    The ``elph_dir`` folder should only be copied or symlinked if ``electron_phonon`` is set in the ``INPUTPH``
    namelist.
    """
    entry_point_name = 'quantumespresso.ph'

    inputs = generate_inputs_ph()
    inputs['parent_folder'] = generate_remote_data(fixture_localhost, str(tmp_path), entry_point_name)
    inputs['settings'] = orm.Dict({'parent_folder_symlink': symlink})

    parameters = {'INPUTPH': {}}
    if electron_phonon is not None:
        parameters['INPUTPH']['electron_phonon'] = electron_phonon
    inputs['parameters'] = orm.Dict(parameters)

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    remote_list = calc_info.remote_symlink_list if symlink else calc_info.remote_copy_list
    targets = [target for _, _, target in remote_list]

    assert (PhCalculation._FOLDER_ELECTRON_PHONON in targets) is (electron_phonon is not None)


def test_serialize_builder(generate_inputs_ph, data_regression, serialize_builder):
    """Test the ``serialize_builder`` fixture using a process builder for the ``PhCalculation``."""
    builder = PhCalculation.get_builder()
    builder._update(**generate_inputs_ph())
    data_regression.check(serialize_builder(builder))


def test_parameters_validation():
    """Test the validation of the `parameters` input."""
    import pytest

    builder = PhCalculation.get_builder()

    parameters = {'inputph': {'tr2_ph': 1.0e-8}}

    with pytest.warns(UserWarning, match="'inputph' should be UPPERCASE"):
        builder.parameters = parameters

    assert builder.parameters.get_dict() == {'INPUTPH': {'tr2_ph': 1.0e-8}}

    with pytest.raises(ValueError, match="'inputph' should be UPPERCASE"):
        builder.parameters = orm.Dict(parameters).store()
