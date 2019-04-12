# pylint: disable=redefined-outer-name
"""Initialise a text database and profile for pytest."""
from __future__ import absolute_import

import shutil
import tempfile

import pytest

from aiida.manage.fixtures import fixture_manager


@pytest.fixture(scope='session')
def fixture_environment():
    """Setup a complete AiiDA test environment, with configuration, profile, database and repository."""
    with fixture_manager() as manager:
        yield manager


@pytest.fixture(scope='session')
def fixture_work_directory():
    """Return a temporary folder that can be used as for example a computer's work directory."""
    dirpath = tempfile.mkdtemp()
    yield dirpath
    shutil.rmtree(dirpath)


@pytest.fixture(scope='function')
def fixture_computer_localhost(fixture_work_directory):
    """Return a `Computer` instance mocking a localhost setup."""
    from aiida.orm import Computer
    computer = Computer(
        name='localhost',
        hostname='localhost',
        transport_type='local',
        scheduler_type='direct',
        workdir=fixture_work_directory).store()
    yield computer


@pytest.fixture(scope='function')
def fixture_database(fixture_environment):
    """Clear the database after each test."""
    yield
    fixture_environment.reset_db()


@pytest.fixture
def generate_calc_job_node():
    """Fixture to generate a mock `CalcJobNode` for testing parsers."""

    def _generate_calc_job_node(entry_point_name, computer, test_name, inputs=None, attributes=None):
        """Fixture to generate a mock `CalcJobNode` for testing parsers.

        :param entry_point_name: entry point name of the calculation class
        :param computer: a `Computer` instance
        :param test_name: relative path of directory with test output files in the `fixtures/{entry_point_name}` folder
        :param attributes: any optional attributes to set on the node
        :return: `CalcJobNode` instance with an attached `FolderData` as the `retrieved` node
        """
        import os
        from aiida import orm
        from aiida.common.links import LinkType
        from aiida.plugins.entry_point import format_entry_point_string

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        node = orm.CalcJobNode(computer=computer, process_type=entry_point)
        node.set_attribute('input_filename', 'aiida.in')
        node.set_attribute('output_filename', 'aiida.out')
        node.set_attribute('error_filename', 'aiida.err')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if attributes:
            node.set_attributes(attributes)

        node.store()

        basepath = os.path.dirname(os.path.abspath(__file__))
        filepath = os.path.join(basepath, 'parsers', 'fixtures', entry_point_name[len('quantumespresso.'):], test_name)

        retrieved = orm.FolderData()
        retrieved.put_object_from_tree(filepath)
        retrieved.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        if inputs:
            for link_label, input_node in inputs.items():
                node.add_incoming(input_node, link_type=LinkType.INPUT_CALC, link_label=link_label)

        return node

    return _generate_calc_job_node


@pytest.fixture
def generate_parser():
    """Fixture to load a parser class for testing parsers."""

    def _generate_parser(entry_point_name):
        """Fixture to load a parser class for testing parsers.

        :param entry_point_name: entry point name of the parser class
        :return: the `Parser` sub class
        """
        from aiida.plugins import ParserFactory
        return ParserFactory(entry_point_name)

    return _generate_parser


@pytest.fixture
def generate_inputs():
    """Fixture to define some basic inputs nodes that are used."""

    def _generate_inputs(entry_point_name):
        """Fixture to load a parser class for testing parsers.

        :param entry_point_name: entry point name of the parser class
        :return: the `Parser` sub class
        """
        from aiida import orm
        from aiida.common import AttributeDict

        if entry_point_name == 'quantumespresso.pw':
            structure = orm.StructureData()
            parameters = {
                'CONTROL': {
                    'calculation': 'scf'
                },
                'SYSTEM': {
                    'ecutrho': 240.0,
                    'ecutwfc': 30.0
                }
            }
            kpoints = orm.KpointsData()
            kpoints.set_cell_from_structure(structure)
            kpoints.set_kpoints_mesh_from_density(0.15)

            inputs = AttributeDict({
                'structure': structure,
                'kpoints': kpoints,
                'parameters': orm.Dict(dict=parameters),
                'settings': orm.Dict()
            })

        return inputs

    return _generate_inputs
