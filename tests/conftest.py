# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name,too-many-statements,too-many-lines
"""Initialise a text database and profile for pytest."""
from collections.abc import Mapping
import io
import os
import pathlib
from pathlib import Path
import shutil
import tempfile

import pytest

pytest_plugins = ['aiida.manage.tests.pytest_fixtures']  # pylint: disable=invalid-name


@pytest.fixture(scope='session')
def filepath_tests():
    """Return the absolute filepath of the `tests` folder.

    .. warning:: if this file moves with respect to the `tests` folder, the implementation should change.

    :return: absolute filepath of `tests` folder which is the basepath for all test resources.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def filepath_fixtures(filepath_tests):
    """Return the absolute filepath to the directory containing the file `fixtures`."""
    return os.path.join(filepath_tests, 'fixtures')


@pytest.fixture(scope='function')
def fixture_sandbox():
    """Return a `SandboxFolder`."""
    from aiida.common.folders import SandboxFolder
    with SandboxFolder() as folder:
        yield folder


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost


@pytest.fixture
def fixture_code(fixture_localhost):
    """Return an ``InstalledCode`` instance configured to run calculations of given entry point on localhost."""

    def _fixture_code(entry_point_name):
        from aiida.common import exceptions
        from aiida.orm import InstalledCode, load_code

        label = f'test.{entry_point_name}'

        try:
            return load_code(label=label)
        except exceptions.NotExistent:
            return InstalledCode(
                label=label,
                computer=fixture_localhost,
                filepath_executable='/bin/true',
                default_calc_job_plugin=entry_point_name,
            )

    return _fixture_code


@pytest.fixture
def serialize_builder():
    """Serialize the given process builder into a dictionary with nodes turned into their value representation.

    :param builder: the process builder to serialize
    :return: dictionary
    """

    def serialize_data(data):
        # pylint: disable=too-many-return-statements
        from aiida.orm import AbstractCode, BaseType, Data, Dict, KpointsData, List, RemoteData, SinglefileData
        from aiida.plugins import DataFactory

        StructureData = DataFactory('core.structure')
        UpfData = DataFactory('pseudo.upf')

        if isinstance(data, dict):
            return {key: serialize_data(value) for key, value in data.items()}

        if isinstance(data, BaseType):
            return data.value

        if isinstance(data, AbstractCode):
            return data.full_label

        if isinstance(data, Dict):
            return data.get_dict()

        if isinstance(data, List):
            return data.get_list()

        if isinstance(data, StructureData):
            return data.get_formula()

        if isinstance(data, UpfData):
            return f'{data.element}<md5={data.md5}>'

        if isinstance(data, RemoteData):
            # For `RemoteData` we compute the hash of the repository. The value returned by `Node._get_hash` is not
            # useful since it includes the hash of the absolute filepath and the computer UUID which vary between tests
            return data.base.repository.hash()

        if isinstance(data, KpointsData):
            try:
                return data.get_kpoints()
            except AttributeError:
                return data.get_kpoints_mesh()

        if isinstance(data, SinglefileData):
            return data.get_content()

        if isinstance(data, Data):
            return data.base.caching._get_hash()  # pylint: disable=protected-access

        return data

    def _serialize_builder(builder):
        return serialize_data(builder._inputs(prune=True))  # pylint: disable=protected-access

    return _serialize_builder


@pytest.fixture(scope='session', autouse=True)
def sssp(aiida_profile, generate_upf_data):
    """Create the SSSP pseudo potential families from scratch."""
    from aiida.common.constants import elements
    from aiida.plugins import GroupFactory

    aiida_profile.clear_profile()

    SsspFamily = GroupFactory('pseudo.family.sssp')

    cutoffs = {}
    stringency = 'standard'

    for label, cutoff_values in zip(('SSSP/1.3/PBEsol/precision', 'SSSP/1.3/PBEsol/efficiency'),
                                    ((40.0, 320.0), (30.0, 240.0))):
        with tempfile.TemporaryDirectory() as dirpath:
            for values in elements.values():

                element = values['symbol']

                actinides = ('Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')

                if element in actinides:
                    continue

                upf = generate_upf_data(element)
                dirpath = pathlib.Path(dirpath)
                filename = dirpath / f'{element}.upf'

                with open(filename, 'w+b') as handle:
                    with upf.open(mode='rb') as source:
                        handle.write(source.read())
                        handle.flush()

                cutoffs[element] = {
                    'cutoff_wfc': cutoff_values[0],
                    'cutoff_rho': cutoff_values[1],
                }

            family = SsspFamily.create_from_folder(dirpath, label)

        family.set_cutoffs(cutoffs, stringency, unit='Ry')

    return family


@pytest.fixture
def generate_calc_job():
    """Fixture to construct a new `CalcJob` instance and call `prepare_for_submission` for testing `CalcJob` classes.

    The fixture will return the `CalcInfo` returned by `prepare_for_submission` and the temporary folder that was passed
    to it, into which the raw input files will have been written.
    """

    def _generate_calc_job(folder, entry_point_name, inputs=None):
        """Fixture to generate a mock `CalcInfo` for testing calculation jobs."""
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager
        from aiida.plugins import CalculationFactory

        manager = get_manager()
        runner = manager.get_runner()

        process_class = CalculationFactory(entry_point_name)
        process = instantiate_process(runner, process_class, **inputs)

        calc_info = process.prepare_for_submission(folder)

        return calc_info

    return _generate_calc_job


@pytest.fixture
def generate_calc_job_node(fixture_localhost):
    """Fixture to generate a mock `CalcJobNode` for testing parsers."""

    def flatten_inputs(inputs, prefix=''):
        """Flatten inputs recursively like :meth:`aiida.engine.processes.process::Process._flatten_inputs`."""
        flat_inputs = []
        for key, value in inputs.items():
            if isinstance(value, Mapping):
                flat_inputs.extend(flatten_inputs(value, prefix=prefix + key + '__'))
            else:
                flat_inputs.append((prefix + key, value))
        return flat_inputs

    def _generate_calc_job_node(
        entry_point_name='base', computer=None, test_name=None, inputs=None, attributes=None, retrieve_temporary=None
    ):
        """Fixture to generate a mock `CalcJobNode` for testing parsers.

        :param entry_point_name: entry point name of the calculation class
        :param computer: a `Computer` instance
        :param test_name: relative path of directory with test output files in the `fixtures/{entry_point_name}` folder.
        :param inputs: any optional nodes to add as input links to the corrent CalcJobNode
        :param attributes: any optional attributes to set on the node
        :param retrieve_temporary: optional tuple of an absolute filepath of a temporary directory and a list of
            filenames that should be written to this directory, which will serve as the `retrieved_temporary_folder`.
            For now this only works with top-level files and does not support files nested in directories.
        :return: `CalcJobNode` instance with an attached `FolderData` as the `retrieved` node.
        """
        from aiida import orm
        from aiida.common import LinkType
        from aiida.plugins.entry_point import format_entry_point_string

        if computer is None:
            computer = fixture_localhost

        filepath_folder = None

        if test_name is not None:
            basepath = os.path.dirname(os.path.abspath(__file__))
            filename = os.path.join(entry_point_name[len('quantumespresso.'):], test_name)
            filepath_folder = os.path.join(basepath, 'parsers', 'fixtures', filename)
            filepath_input = os.path.join(filepath_folder, 'aiida.in')

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        node = orm.CalcJobNode(computer=computer, process_type=entry_point)
        node.base.attributes.set('input_filename', 'aiida.in')
        node.base.attributes.set('output_filename', 'aiida.out')
        node.base.attributes.set('error_filename', 'aiida.err')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if attributes:
            node.base.attributes.set_many(attributes)

        if filepath_folder:
            from qe_tools.exceptions import ParsingError

            from aiida_quantumespresso.tools.pwinputparser import PwInputFile
            try:
                with open(filepath_input, 'r', encoding='utf-8') as input_file:
                    parsed_input = PwInputFile(input_file.read())
            except (ParsingError, FileNotFoundError):
                pass
            else:
                inputs['structure'] = parsed_input.get_structuredata()
                inputs['parameters'] = orm.Dict(parsed_input.namelists)

        if inputs:
            metadata = inputs.pop('metadata', {})
            options = metadata.get('options', {})

            for name, option in options.items():
                node.set_option(name, option)

            for link_label, input_node in flatten_inputs(inputs):
                input_node.store()
                node.base.links.add_incoming(input_node, link_type=LinkType.INPUT_CALC, link_label=link_label)

        node.store()

        if retrieve_temporary:
            dirpath, filenames = retrieve_temporary
            dirpath = Path(dirpath)
            filepaths = []
            for filename in filenames:
                filepaths.extend(Path(filepath_folder).glob(filename))

            for filepath in filepaths:
                shutil.copy(filepath, dirpath / filepath.name)

        if filepath_folder:
            retrieved = orm.FolderData()
            retrieved.base.repository.put_object_from_tree(filepath_folder)

            # Remove files that are supposed to be only present in the retrieved temporary folder
            if retrieve_temporary:
                for filepath in filepaths:
                    retrieved.delete_object(filepath.name)

            retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
            retrieved.store()

            remote_folder = orm.RemoteData(computer=computer, remote_path='/tmp')
            remote_folder.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='remote_folder')
            remote_folder.store()

        return node

    return _generate_calc_job_node


@pytest.fixture(scope='session')
def generate_upf_data():
    """Return a `UpfData` instance for the given element a file for which should exist in `tests/fixtures/pseudos`."""

    def _generate_upf_data(element):
        """Return `UpfData` node."""
        from aiida_pseudo.data.pseudo import UpfData
        content = f'<UPF version="2.0.1"><PP_HEADER\nelement="{element}"\nz_valence="4.0"\n/></UPF>\n'
        stream = io.BytesIO(content.encode('utf-8'))
        return UpfData(stream, filename=f'{element}.upf')

    return _generate_upf_data


@pytest.fixture(scope='session')
def generate_xy_data():
    """Return an ``XyData`` instance."""

    def _generate_xy_data():
        """Return an ``XyData`` node."""
        from aiida.orm import XyData
        import numpy as np

        xvals = [1, 2, 3]
        yvals = [10, 20, 30]
        xlabel = 'X'
        ylabel = 'Y'
        xunits = 'n/a'
        yunits = 'n/a'

        xy_node = XyData()
        xy_node.set_x(np.array(xvals), xlabel, xunits)
        xy_node.set_y([np.array(yvals)], [ylabel], [yunits])

        return xy_node

    return _generate_xy_data


@pytest.fixture
def generate_structure():
    """Return a ``StructureData`` representing either bulk silicon or a water molecule."""

    def _generate_structure(structure_id='silicon'):
        """Return a ``StructureData`` representing bulk silicon or a snapshot of a single water molecule dynamics.

        :param structure_id: identifies the ``StructureData`` you want to generate. Either 'silicon' or 'water'.
        """
        from aiida.orm import StructureData

        if structure_id.startswith('silicon'):
            name1 = 'Si0' if structure_id.endswith('kinds') else 'Si'
            name2 = 'Si1' if structure_id.endswith('kinds') else 'Si'
            param = 5.43
            cell = [[param / 2., param / 2., 0], [param / 2., 0, param / 2.], [0, param / 2., param / 2.]]
            structure = StructureData(cell=cell)
            structure.append_atom(position=(0., 0., 0.), symbols='Si', name=name1)
            structure.append_atom(position=(param / 4., param / 4., param / 4.), symbols='Si', name=name2)
        elif structure_id == 'cobalt-prim':
            cell = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]
            structure = StructureData(cell=cell)
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols='Co', name='Co')
        elif structure_id == 'water':
            structure = StructureData(cell=[[5.29177209, 0., 0.], [0., 5.29177209, 0.], [0., 0., 5.29177209]])
            structure.append_atom(position=[12.73464656, 16.7741411, 24.35076238], symbols='H', name='H')
            structure.append_atom(position=[-29.3865565, 9.51707929, -4.02515904], symbols='H', name='H')
            structure.append_atom(position=[1.04074437, -1.64320127, -1.27035021], symbols='O', name='O')
        elif structure_id == 'uranium':
            param = 5.43
            cell = [[param / 2., param / 2., 0], [param / 2., 0, param / 2.], [0, param / 2., param / 2.]]
            structure = StructureData(cell=cell)
            structure.append_atom(position=(0., 0., 0.), symbols='U', name='U')
            structure.append_atom(position=(param / 4., param / 4., param / 4.), symbols='U', name='U')
        elif structure_id == '2D-xy-arsenic':
            cell = [[3.61, 0, 0], [-1.80, 3.13, 0], [0, 0, 21.3]]
            structure = StructureData(cell=cell, pbc=(True, True, False))
            structure.append_atom(position=(1.804, 1.042, 11.352), symbols='As', name='As')
            structure.append_atom(position=(0, 2.083, 9.960), symbols='As', name='As')
        elif structure_id == '1D-x-carbon':
            cell = [[4.2, 0, 0], [0, 20, 0], [0, 0, 20]]
            structure = StructureData(cell=cell, pbc=(True, False, False))
            structure.append_atom(position=(0, 0, 0), symbols='C', name='C')
        elif structure_id == '1D-y-carbon':
            cell = [[20, 0, 0], [0, 4.2, 0], [0, 0, 20]]
            structure = StructureData(cell=cell, pbc=(False, True, False))
            structure.append_atom(position=(0, 0, 0), symbols='C', name='C')
        elif structure_id == '1D-z-carbon':
            cell = [[20, 0, 0], [0, 20, 0], [0, 0, 4.2]]
            structure = StructureData(cell=cell, pbc=(False, False, True))
            structure.append_atom(position=(0, 0, 0), symbols='C', name='C')
        else:
            raise KeyError(f'Unknown structure_id="{structure_id}"')
        return structure

    return _generate_structure


@pytest.fixture
def generate_structure_from_kinds():
    """Return a dummy `StructureData` instance with the specified kind names."""

    def _generate_structure_from_kinds(site_kind_names):
        """Return a dummy `StructureData` instance with the specified kind names."""
        import re

        from aiida import orm

        if not isinstance(site_kind_names, (list, tuple)):
            site_kind_names = (site_kind_names,)

        structure = orm.StructureData(cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        for kind_name in site_kind_names:
            structure.append_atom(name=kind_name, symbols=re.sub('[0-9]', '', kind_name), position=(0., 0., 0.))

        return structure

    return _generate_structure_from_kinds


@pytest.fixture
def generate_kpoints_mesh():
    """Return a `KpointsData` node."""

    def _generate_kpoints_mesh(npoints):
        """Return a `KpointsData` with a mesh of npoints in each direction."""
        from aiida.orm import KpointsData

        kpoints = KpointsData()
        kpoints.set_kpoints_mesh([npoints] * 3)

        return kpoints

    return _generate_kpoints_mesh


@pytest.fixture(scope='session')
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
def generate_remote_data():
    """Return a `RemoteData` node."""

    def _generate_remote_data(computer, remote_path, entry_point_name=None):
        """Return a `RemoteData` node."""
        from aiida.common.links import LinkType
        from aiida.orm import CalcJobNode, RemoteData
        from aiida.plugins.entry_point import format_entry_point_string

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        remote = RemoteData(remote_path=remote_path)
        remote.computer = computer

        if entry_point_name is not None:
            creator = CalcJobNode(computer=computer, process_type=entry_point)
            creator.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
            remote.base.links.add_incoming(creator, link_type=LinkType.CREATE, link_label='remote_folder')
            creator.store()

        return remote

    return _generate_remote_data


@pytest.fixture
def generate_bands_data():
    """Return a `BandsData` node."""

    def _generate_bands_data():
        """Return a `BandsData` instance with some basic `kpoints` and `bands` arrays."""
        from aiida.plugins import DataFactory
        import numpy
        BandsData = DataFactory('core.array.bands')  #pylint: disable=invalid-name
        bands_data = BandsData()

        bands_data.set_kpoints(numpy.array([[0., 0., 0.], [0.625, 0.25, 0.625]]))

        bands_data.set_bands(
            numpy.array([[-5.64024889, 6.66929678, 6.66929678, 6.66929678, 8.91047649],
                         [-1.71354964, -0.74425095, 1.82242466, 3.98697455, 7.37979746]]),
            units='eV'
        )

        return bands_data

    return _generate_bands_data


@pytest.fixture
def generate_workchain():
    """Generate an instance of a `WorkChain`."""

    def _generate_workchain(entry_point, inputs):
        """Generate an instance of a `WorkChain` with the given entry point and inputs.

        :param entry_point: entry point name of the work chain subclass.
        :param inputs: inputs to be passed to process construction.
        :return: a `WorkChain` instance.
        """
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager
        from aiida.plugins import WorkflowFactory

        process_class = WorkflowFactory(entry_point)
        runner = get_manager().get_runner()
        process = instantiate_process(runner, process_class, **inputs)

        return process

    return _generate_workchain


@pytest.fixture
def generate_force_constants_data(filepath_tests):
    """Generate a ``ForceConstantsData`` node."""
    from aiida_quantumespresso.data.force_constants import ForceConstantsData
    filepath = os.path.join(filepath_tests, 'calculations', 'fixtures', 'matdyn', 'default', 'force_constants.dat')
    return ForceConstantsData(filepath)


@pytest.fixture
def generate_inputs_matdyn(fixture_code, generate_kpoints_mesh, generate_force_constants_data):
    """Generate default inputs for a `MatdynCalculation."""

    def _generate_inputs_matdyn():
        """Generate default inputs for a `MatdynCalculation."""
        from aiida_quantumespresso.utils.resources import get_default_options

        inputs = {
            'code': fixture_code('quantumespresso.matdyn'),
            'force_constants': generate_force_constants_data,
            'kpoints': generate_kpoints_mesh(2),
            'metadata': {
                'options': get_default_options()
            }
        }

        return inputs

    return _generate_inputs_matdyn


@pytest.fixture
def generate_inputs_q2r(fixture_sandbox, fixture_localhost, fixture_code, generate_remote_data):
    """Generate default inputs for a `Q2rCalculation."""

    def _generate_inputs_q2r():
        """Generate default inputs for a `Q2rCalculation."""
        from aiida_quantumespresso.utils.resources import get_default_options

        inputs = {
            'code': fixture_code('quantumespresso.q2r'),
            'parent_folder': generate_remote_data(fixture_localhost, fixture_sandbox.abspath, 'quantumespresso.ph'),
            'metadata': {
                'options': get_default_options()
            }
        }

        return inputs

    return _generate_inputs_q2r


@pytest.fixture
def generate_inputs_bands(fixture_sandbox, fixture_localhost, fixture_code, generate_remote_data):
    """Generate default inputs for a `BandsCalculation."""

    def _generate_inputs_bands():
        """Generate default inputs for a `BandsCalculation."""
        from aiida_quantumespresso.utils.resources import get_default_options

        inputs = {
            'code': fixture_code('quantumespresso.bands'),
            'parent_folder': generate_remote_data(fixture_localhost, fixture_sandbox.abspath, 'quantumespresso.pw'),
            'metadata': {
                'options': get_default_options()
            }
        }

        return inputs

    return _generate_inputs_bands


@pytest.fixture
def generate_inputs_ph(
    generate_calc_job_node, generate_structure, fixture_localhost, fixture_code, generate_kpoints_mesh
):
    """Generate default inputs for a `PhCalculation."""

    def _generate_inputs_ph(with_output_structure=False):
        """Generate default inputs for a `PhCalculation.

        :param with_output_structure: whether the PwCalculation has a StructureData in its outputs.
            This is needed to test some PhBaseWorkChain logics.
        """
        from aiida.common import LinkType
        from aiida.orm import Dict, RemoteData

        from aiida_quantumespresso.utils.resources import get_default_options

        pw_node = generate_calc_job_node(
            entry_point_name='quantumespresso.pw', inputs={
                'parameters': Dict(),
                'structure': generate_structure()
            }
        )
        remote_folder = RemoteData(computer=fixture_localhost, remote_path='/tmp')
        remote_folder.base.links.add_incoming(pw_node, link_type=LinkType.CREATE, link_label='remote_folder')
        remote_folder.store()
        parent_folder = pw_node.outputs.remote_folder

        if with_output_structure:
            structure = generate_structure()
            structure.base.links.add_incoming(pw_node, link_type=LinkType.CREATE, link_label='output_structure')
            structure.store()

        inputs = {
            'code': fixture_code('quantumespresso.ph'),
            'parent_folder': parent_folder,
            'qpoints': generate_kpoints_mesh(2),
            'parameters': Dict({'INPUTPH': {}}),
            'metadata': {
                'options': get_default_options()
            }
        }

        return inputs

    return _generate_inputs_ph


@pytest.fixture
def generate_inputs_pw(fixture_code, generate_structure, generate_kpoints_mesh, generate_upf_data):
    """Generate default inputs for a `PwCalculation."""

    def _generate_inputs_pw():
        """Generate default inputs for a `PwCalculation."""
        from aiida.orm import Dict

        from aiida_quantumespresso.utils.resources import get_default_options

        parameters = Dict({
            'CONTROL': {
                'calculation': 'scf'
            },
            'SYSTEM': {
                'ecutrho': 240.0,
                'ecutwfc': 30.0
            },
            'ELECTRONS': {
                'electron_maxstep': 60,
            }
        })
        structure = generate_structure()
        inputs = {
            'code': fixture_code('quantumespresso.pw'),
            'structure': generate_structure(),
            'kpoints': generate_kpoints_mesh(2),
            'parameters': parameters,
            'pseudos': {kind: generate_upf_data(kind) for kind in structure.get_kind_names()},
            'metadata': {
                'options': get_default_options()
            }
        }
        return inputs

    return _generate_inputs_pw


@pytest.fixture
def generate_inputs_cp(fixture_code, generate_structure, generate_upf_data):
    """Generate default inputs for a CpCalculation."""

    def _generate_inputs_cp(autopilot=False):
        """Generate default inputs for a CpCalculation."""
        from aiida.orm import Dict

        from aiida_quantumespresso.utils.resources import get_default_options

        inputs = {
            'code': fixture_code('quantumespresso.cp'),
            'structure': generate_structure(),
            'parameters': Dict({
                'CONTROL': {
                    'calculation': 'cp'
                },
                'SYSTEM': {
                    'ecutrho': 240.0,
                    'ecutwfc': 30.0
                }
            }),
            'pseudos': {
                'Si': generate_upf_data('Si')
            },
            'metadata': {
                'options': get_default_options()
            }
        }
        if autopilot:
            inputs['settings'] = Dict({
                'AUTOPILOT': [{
                    'onstep': 2,
                    'what': 'dt',
                    'newvalue': 42.0
                }, {
                    'onstep': 3,
                    'what': 'dt',
                    'newvalue': 42.42
                }]
            })

        return inputs

    return _generate_inputs_cp


@pytest.fixture
def generate_inputs_xspectra(
    fixture_sandbox,
    fixture_localhost,
    fixture_code,
    generate_remote_data,
    generate_kpoints_mesh,
):
    """Generate inputs for an ``XspectraCalculation``."""

    def _generate_inputs_xspectra():
        from aiida.orm import Dict, SinglefileData

        from aiida_quantumespresso.utils.resources import get_default_options

        parameters = {
            'INPUT_XSPECTRA': {
                'calculation': 'xanes_dipole',
            },
        }

        inputs = {
            'code':
            fixture_code('quantumespresso.xspectra'),
            'parameters':
            Dict(parameters),
            'parent_folder':
            generate_remote_data(fixture_localhost, fixture_sandbox.abspath, 'quantumespresso.pw'),
            'core_wfc_data':
            SinglefileData(
                io.StringIO(
                    '# number of core states 3 =  1 0;  2 0;'
                    '\n6.51344e-05 6.615743462459999e-3'
                    '\n6.59537e-05 6.698882211449999e-3'
                )
            ),
            'kpoints':
            generate_kpoints_mesh(2),
            'metadata': {
                'options': get_default_options()
            }
        }

        return inputs

    return _generate_inputs_xspectra


@pytest.fixture
def generate_workchain_pw(generate_workchain, generate_inputs_pw, generate_calc_job_node):
    """Generate an instance of a ``PwBaseWorkChain``."""

    def _generate_workchain_pw(exit_code=None, inputs=None, return_inputs=False, pw_outputs=None):
        """Generate an instance of a ``PwBaseWorkChain``.

        :param exit_code: exit code for the ``PwCalculation``.
        :param inputs: inputs for the ``PwBaseWorkChain``.
        :param return_inputs: return the inputs of the ``PwBaseWorkChain``.
        :param pw_outputs: ``dict`` of outputs for the ``PwCalculation``. The keys must correspond to the link labels
            and the values to the output nodes.
        """
        from aiida.common import LinkType
        from aiida.orm import Dict
        from plumpy import ProcessState

        entry_point = 'quantumespresso.pw.base'

        if inputs is None:
            pw_inputs = generate_inputs_pw()
            kpoints = pw_inputs.pop('kpoints')
            inputs = {'pw': pw_inputs, 'kpoints': kpoints}

        if return_inputs:
            return inputs

        process = generate_workchain(entry_point, inputs)

        pw_node = generate_calc_job_node(inputs={'parameters': Dict()})
        process.ctx.iteration = 1
        process.ctx.children = [pw_node]

        if pw_outputs is not None:
            for link_label, output_node in pw_outputs.items():
                output_node.base.links.add_incoming(pw_node, link_type=LinkType.CREATE, link_label=link_label)
                output_node.store()

        if exit_code is not None:
            pw_node.set_process_state(ProcessState.FINISHED)
            pw_node.set_exit_status(exit_code.status)

        return process

    return _generate_workchain_pw


@pytest.fixture
def generate_workchain_ph(generate_workchain, generate_inputs_ph, generate_calc_job_node):
    """Generate an instance of a `PhBaseWorkChain`."""

    def _generate_workchain_ph(exit_code=None, inputs=None, return_inputs=False):
        from plumpy import ProcessState

        entry_point = 'quantumespresso.ph.base'

        if inputs is None:
            ph_inputs = generate_inputs_ph()
            qpoints = ph_inputs.pop('qpoints')
            inputs = {'ph': ph_inputs, 'qpoints': qpoints}

        if return_inputs:
            return inputs

        process = generate_workchain(entry_point, inputs)

        if exit_code is not None:
            node = generate_calc_job_node()
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.iteration = 1
            process.ctx.children = [node]

        return process

    return _generate_workchain_ph


@pytest.fixture
def generate_workchain_pdos(generate_workchain, generate_inputs_pw, fixture_code):
    """Generate an instance of a `PdosWorkChain`."""

    def _generate_workchain_pdos(emin=None, emax=None, energy_range_vs_fermi=None):
        from aiida.orm import Bool, Dict, List

        from aiida_quantumespresso.utils.resources import get_default_options

        entry_point = 'quantumespresso.pdos'

        scf_pw_inputs = generate_inputs_pw()
        kpoints = scf_pw_inputs.pop('kpoints')
        structure = scf_pw_inputs.pop('structure')
        scf = {'pw': scf_pw_inputs, 'kpoints': kpoints}

        nscf_pw_inputs = generate_inputs_pw()
        nscf_pw_inputs.pop('kpoints')
        nscf_pw_inputs.pop('structure')
        nscf_pw_inputs['parameters']['CONTROL']['calculation'] = 'nscf'
        nscf_pw_inputs['parameters']['SYSTEM']['occupations'] = 'tetrahedra'
        nscf_pw_inputs['parameters']['SYSTEM']['nosym'] = True

        nscf = {'pw': nscf_pw_inputs, 'kpoints': kpoints}

        dos_params = {
            'DOS': {
                'DeltaE': 0.01,
            }
        }
        projwfc_params = {'PROJWFC': {'DeltaE': 0.01, 'ngauss': 0, 'degauss': 0.01}}

        if emin and emax:
            dos_params['DOS'].update({'Emin': emin, 'Emax': emax})
            projwfc_params['PROJWFC'].update({'Emin': emin, 'Emax': emax})

        dos = {
            'code': fixture_code('quantumespresso.dos'),
            'parameters': Dict(dos_params),
            'metadata': {
                'options': get_default_options()
            }
        }
        projwfc = {
            'code': fixture_code('quantumespresso.projwfc'),
            'parameters': Dict(projwfc_params),
            'metadata': {
                'options': get_default_options()
            }
        }
        inputs = {
            'structure': structure,
            'scf': scf,
            'nscf': nscf,
            'dos': dos,
            'projwfc': projwfc,
            'dry_run': Bool(True)
        }
        if energy_range_vs_fermi:
            inputs.update({'energy_range_vs_fermi': List(energy_range_vs_fermi)})

        return generate_workchain(entry_point, inputs)

    return _generate_workchain_pdos


@pytest.fixture
def generate_workchain_xspectra(generate_workchain, generate_inputs_xspectra, generate_calc_job_node):
    """Generate an instance of a `XspectraBaseWorkChain`."""

    def _generate_workchain_xspectra(exit_code=None, inputs=None, return_inputs=False, xspectra_outputs=None):
        from aiida.common import LinkType
        from plumpy import ProcessState

        entry_point = 'quantumespresso.xspectra.base'

        if inputs is None:
            inputs_xspectra = generate_inputs_xspectra()
            kpoints = inputs_xspectra.pop('kpoints')
            inputs = {'xspectra': inputs_xspectra, 'kpoints': kpoints}

        if return_inputs:
            return inputs

        process = generate_workchain(entry_point, inputs)
        xspectra_node = generate_calc_job_node()
        process.ctx.iteration = 1
        process.ctx.children = [xspectra_node]

        if xspectra_outputs is not None:
            for link_label, output_node in xspectra_outputs.items():
                output_node.base.links.add_incoming(xspectra_node, link_type=LinkType.CREATE, link_label=link_label)
                output_node.store()

        if exit_code is not None:
            xspectra_node.set_process_state(ProcessState.FINISHED)
            xspectra_node.set_exit_status(exit_code.status)

        return process

    return _generate_workchain_xspectra


@pytest.fixture
def generate_workchain_xps(generate_inputs_pw, generate_workchain, generate_upf_data):
    """Generate an instance of a `XpsWorkChain`."""

    def _generate_workchain_xps():
        from aiida.orm import Bool, List, Str

        entry_point = 'quantumespresso.xps'

        scf_pw_inputs = generate_inputs_pw()
        kpoints = scf_pw_inputs.pop('kpoints')
        structure = scf_pw_inputs.pop('structure')
        ch_scf = {'pw': scf_pw_inputs, 'kpoints': kpoints}

        inputs = {
            'structure': structure,
            'ch_scf': ch_scf,
            'dry_run': Bool(True),
            'elements_list': List(['Si']),
            'abs_atom_marker': Str('X'),
            'core_hole_pseudos': {
                'Si': generate_upf_data('Si')
            },
            'gipaw_pseudos': {
                'Si': generate_upf_data('Si')
            },
        }

        return generate_workchain(entry_point, inputs)

    return _generate_workchain_xps
