# -*- coding: utf-8 -*-
"""Tests for the `EpwCalculation` class."""
from aiida import orm
from aiida.plugins import CalculationFactory
from aiida.common.links import LinkType
from aiida_quantumespresso.utils.resources import get_default_options

EPWCALC = CalculationFactory('quantumespresso.epw')


def test_epw_default(
    fixture_localhost, fixture_sandbox, generate_calc_job, generate_remote_data, fixture_code, generate_kpoints_mesh,
    file_regression, tmpdir
):
    """Test a default `EpwCalculation`."""
    entry_point_name = 'quantumespresso.epw'

    parameters = {
        'INPUTEPW': {
            'nbndsub': 8,
            'elph': True,  # default is false
            'epbwrite': True,  # default is false
            'epwwrite': True,  # default is false
            'proj(1)': 'Si : sp3',
            'elecselfen': True,
            'wannierize': True,
            'dvscf_dir': './save/',
            'dis_win_max': 18,
            'dis_froz_max': 8.5
        }
    }

    parent_pw = generate_remote_data(
        fixture_localhost,
        str(tmpdir),
        'quantumespresso.pw',
    )

    parent_ph = generate_remote_data(
        fixture_localhost,
        str(tmpdir),
        'quantumespresso.ph',
    )

    parent_ph.set_remote_path('dummy_path')

    # Dummy q-point
    qibz = {
        'dynamical_matrix_1': {
            'q_point': [0.0, 0.0, 0.0],
        },
        'dynamical_matrix_2': {
            'q_point': [0.0, 0.5, 0.5],
        },
        'dynamical_matrix_3': {
            'q_point': [0.5, 0.5, 0.5],
        }
    }

    parameters2 = orm.Dict(dict=qibz)
    parameters2.add_incoming(parent_ph.creator, link_label='output_parameters', link_type=LinkType.CREATE)
    parameters2.store()

    inputs = {
        'code': fixture_code(entry_point_name),
        'qpoints': generate_kpoints_mesh(2),
        'kpoints': generate_kpoints_mesh(2),
        'qfpoints': generate_kpoints_mesh(2),
        'kfpoints': generate_kpoints_mesh(2),
        'parent_folder_nscf': parent_pw,
        'parent_folder_ph': parent_ph,
        'parameters': orm.Dict(dict=parameters),
        'metadata': {
            'options': get_default_options()
        }
    }

    generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.in')
