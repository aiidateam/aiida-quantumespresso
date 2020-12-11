# -*- coding: utf-8 -*-
import re
import fnmatch

import numpy as np

from aiida.common import LinkType
from aiida.orm import Dict, ProjectionData, BandsData, XyData, CalcJobNode
from aiida.plugins import OrbitalFactory

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw.base import parse_output_base
from .base import Parser


def find_orbitals_from_statelines(out_info_dict):
    """This function reads in all the state_lines, that is, the lines describing which atomic states, taken from the
    pseudopotential, are used for the projection. Then it converts these state_lines into a set of orbitals.

    :param out_info_dict: contains various technical internals useful in parsing
    :return: orbitals, a list of orbitals suitable for setting ProjectionData
    """

    # Format of statelines
    # From PP/src/projwfc.f90: (since Oct. 8 2019)
    #
    # 1000 FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2," (l=",i1)
    # IF (lspinorb) THEN
    # 1001 FORMAT (" j=",f3.1," m_j=",f4.1,")")
    # ELSE IF (noncolin) THEN
    # 1002 FORMAT (" m=",i2," s_z=",f4.1,")")
    # ELSE
    # 1003 FORMAT (" m=",i2,")")
    # ENDIF
    #
    # Before:
    # IF (lspinorb) THEN
    # ...
    # 1000    FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
    #               " (j=",f3.1," l=",i1," m_j=",f4.1,")")
    # ELSE
    # ...
    # 1500    FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
    #               " (l=",i1," m=",i2," s_z=",f4.1,")")
    # ENDIF

    out_file = out_info_dict['out_file']
    atomnum_re = re.compile(r'atom\s*([0-9]+?)[^0-9]')
    element_re = re.compile(r'atom\s*[0-9]+\s*\(\s*([A-Za-z0-9-_]+?)\s*\)')
    if out_info_dict['spinorbit']:
        # spinorbit
        lnum_re = re.compile(r'l=\s*([0-9]+?)[^0-9]')
        jnum_re = re.compile(r'j=\s*([0-9.]+?)[^0-9.]')
        mjnum_re = re.compile(r'm_j=\s*([-0-9.]+?)[^-0-9.]')
    elif not out_info_dict['collinear']:
        # non-collinear
        lnum_re = re.compile(r'l=\s*([0-9]+?)[^0-9]')
        mnum_re = re.compile(r'm=\s*([-0-9]+?)[^-0-9]')
        sznum_re = re.compile(r's_z=\s*([-0-9.]*?)[^-0-9.]')
    else:
        # collinear / no spin
        lnum_re = re.compile(r'l=\s*([0-9]+?)[^0-9]')
        mnum_re = re.compile(r'm=\s*([-0-9]+?)[^-0-9]')
    wfc_lines = out_info_dict['wfc_lines']
    state_lines = [out_file[wfc_line] for wfc_line in wfc_lines]
    state_dicts = []
    for state_line in state_lines:
        try:
            state_dict = {}
            state_dict['atomnum'] = int(atomnum_re.findall(state_line)[0])
            state_dict['atomnum'] -= 1  # to keep with orbital indexing
            state_dict['kind_name'] = element_re.findall(state_line)[0].strip()
            state_dict['angular_momentum'] = int(lnum_re.findall(state_line)[0])
            if out_info_dict['spinorbit']:
                state_dict['total_angular_momentum'] = float(jnum_re.findall(state_line)[0])
                state_dict['magnetic_number'] = float(mjnum_re.findall(state_line)[0])
            else:
                if not out_info_dict['collinear']:
                    state_dict['spin'] = float(sznum_re.findall(state_line)[0])
                state_dict['magnetic_number'] = int(mnum_re.findall(state_line)[0])
                state_dict['magnetic_number'] -= 1  # to keep with orbital indexing
        except ValueError:
            raise QEOutputParsingError('State lines are not formatted in a standard way.')
        state_dicts.append(state_dict)

    # here is some logic to figure out the value of radial_nodes to use
    new_state_dicts = []
    for i in range(len(state_dicts)):
        radial_nodes = 0
        state_dict = state_dicts[i].copy()
        for j in range(i - 1, -1, -1):
            if state_dict == state_dicts[j]:
                radial_nodes += 1
        state_dict['radial_nodes'] = radial_nodes
        new_state_dicts.append(state_dict)
    state_dicts = new_state_dicts

    # here is some logic to assign positions based on the atom_index
    structure = out_info_dict['structure']
    for state_dict in state_dicts:
        site_index = state_dict.pop('atomnum')
        state_dict['position'] = structure.sites[site_index].position

    # here we set the resulting state_dicts to a new set of orbitals
    orbitals = []
    if out_info_dict['spinorbit']:
        OrbitalCls = OrbitalFactory('spinorbithydrogen')
    elif not out_info_dict['collinear']:
        OrbitalCls = OrbitalFactory('noncollinearhydrogen')
    else:
        OrbitalCls = OrbitalFactory('realhydrogen')
    for state_dict in state_dicts:
        orbitals.append(OrbitalCls(**state_dict))

    return orbitals


def spin_dependent_subparser(out_info_dict):
    """This find the projection and bands arrays from the out_file and out_info_dict. Used to handle the different
    possible spin-cases in a convenient manner.

    :param out_info_dict: contains various technical internals useful in parsing
    :return: ProjectionData, BandsData parsed from out_file
    """
    out_file = out_info_dict['out_file']
    spin_down = out_info_dict['spin_down']
    od = out_info_dict  # using a shorter name for convenience
    #   regular expressions needed for later parsing
    WaveFraction1_re = re.compile(r'\=(.*?)\*')  # state composition 1
    WaveFractionremain_re = re.compile(r'\+(.*?)\*')  # state comp 2
    FunctionId_re = re.compile(r'\#(.*?)\]')  # state identity
    # primes arrays for the later parsing
    num_wfc = len(od['wfc_lines'])
    bands = np.zeros([od['k_states'], od['num_bands']])
    projection_arrays = np.zeros([od['k_states'], od['num_bands'], num_wfc])

    try:
        for i in range(od['k_states']):
            if spin_down:
                i += od['k_states']
            # grabs band energy
            for j in range(i * od['num_bands'], (i + 1) * od['num_bands'], 1):
                out_ind = od['e_lines'][j]
                try:
                    # post ~6.3 <6.5 output format "e ="
                    val = out_file[out_ind].split()[2]
                    float(val)
                except ValueError:
                    # pre ~6.3 and 6.5+ output format "==== e("
                    val = out_file[out_ind].split(' = ')[1].split()[0]
                bands[i % od['k_states']][j % od['num_bands']] = val
                #subloop grabs pdos
                wave_fraction = []
                wave_id = []
                for k in range(od['e_lines'][j] + 1, od['psi_lines'][j], 1):
                    out_line = out_file[k]
                    wave_fraction += WaveFraction1_re.findall(out_line)
                    wave_fraction += WaveFractionremain_re.findall(out_line)
                    wave_id += FunctionId_re.findall(out_line)
                if len(wave_id) != len(wave_fraction):
                    raise IndexError
                for l in range(len(wave_id)):
                    wave_id[l] = int(wave_id[l])
                    wave_fraction[l] = float(wave_fraction[l])
                    #sets relevant values in pdos_array
                    projection_arrays[i % od['k_states']][j % od['num_bands']][wave_id[l] - 1] = wave_fraction[l]
    except IndexError:
        raise QEOutputParsingError('the standard out file does not comply with the official documentation.')

    bands_data = BandsData()
    # Attempts to retrieve the kpoints from the parent calc
    parent_calc = out_info_dict['parent_calc']
    try:
        parent_kpoints = parent_calc.get_incoming(link_label_filter='kpoints').one().node
    except ValueError:
        raise QEOutputParsingError('The parent had no input kpoints! Cannot parse from this!')
    try:
        if len(od['k_vect']) != len(parent_kpoints.get_kpoints()):
            raise AttributeError
        bands_data.set_kpointsdata(parent_kpoints)
    except AttributeError:
        bands_data.set_kpoints(od['k_vect'].astype(float))

    bands_data.set_bands(bands, units='eV')

    orbitals = out_info_dict['orbitals']
    if len(orbitals) != np.shape(projection_arrays[0, 0, :])[0]:
        raise QEOutputParsingError(
            'There was an internal parsing error, '
            ' the projection array shape does not agree'
            ' with the number of orbitals'
        )
    projection_data = ProjectionData()
    projection_data.set_reference_bandsdata(bands_data)
    projections = [projection_arrays[:, :, i] for i in range(len(orbitals))]

    # Do the bands_check manually here
    for projection in projections:
        if np.shape(projection) != np.shape(bands):
            raise AttributeError('Projections not the same shape as the bands')

    #insert here some logic to assign pdos to the orbitals
    pdos_arrays = spin_dependent_pdos_subparser(out_info_dict)
    energy_arrays = [out_info_dict['energy']] * len(orbitals)
    projection_data.set_projectiondata(
        orbitals,
        list_of_projections=projections,
        list_of_energy=energy_arrays,
        list_of_pdos=pdos_arrays,
        bands_check=False
    )
    # pdos=pdos_arrays
    return bands_data, projection_data


def natural_sort_key(sort_key, _nsre=re.compile('([0-9]+)')):
    """Pass to ``key`` for ``str.sort`` to achieve natural sorting. For example, ``["2", "11", "1"]`` will be sorted to
    ``["1", "2", "11"]`` instead of ``["1", "11", "2"]``

    :param sort_key: Original key to be processed
    :return: A list of string and integers.
    """
    keys = []
    for text in _nsre.split(sort_key):
        if text.isdigit():
            keys.append(int(text))
        else:
            keys.append(text)
    return keys


def spin_dependent_pdos_subparser(out_info_dict):
    """Finds and labels the pdos arrays associated with the out_info_dict.

    :param out_info_dict: contains various technical internals useful in parsing
    :return: (pdos_name, pdos_array) tuples for all the specific pdos
    """
    spin = out_info_dict['spin']
    collinear = out_info_dict.get('collinear', True)
    spinorbit = out_info_dict.get('spinorbit', False)
    spin_down = out_info_dict['spin_down']
    pdos_atm_array_dict = out_info_dict['pdos_atm_array_dict']
    if spin:
        mult_factor = 2
        if spin_down:
            first_array = 4
        else:
            first_array = 3
    else:
        mult_factor = 1
        first_array = 2
    mf = mult_factor
    fa = first_array
    pdos_file_names = [k for k in pdos_atm_array_dict]
    pdos_file_names.sort(key=natural_sort_key)
    out_arrays = []
    # we can keep the pdos in synch with the projections by relying on the fact
    # both are produced in the same order (thus the sorted file_names)
    for name in pdos_file_names:
        this_array = pdos_atm_array_dict[name]
        if not collinear and not spinorbit:
            # In the non-collinear, non-spinorbit case, the "up"-spin orbitals
            # come first, followed by all "down" orbitals
            for i in range(3, np.shape(this_array)[1], 2):
                out_arrays.append(this_array[:, i])
            for i in range(4, np.shape(this_array)[1], 2):
                out_arrays.append(this_array[:, i])
        else:
            for i in range(fa, np.shape(this_array)[1], mf):
                out_arrays.append(this_array[:, i])

    return out_arrays


class ProjwfcParser(Parser):
    """This class is the implementation of the Parser class for projwfc.x in Quantum Espresso.

    Parses projection arrays that map the projection onto each point in the bands structure, as well as pdos arrays,
    which map the projected density of states onto an energy axis.
    """

    def parse(self, **kwargs):
        """Parses the datafolder, stores results.

        Retrieves projwfc output, and some basic information from the out_file, such as warnings and wall_time
        """
        # Check that the retrieved folder is there
        retrieved = self.retrieved

        # Read standard out
        try:
            filename_stdout = self.node.get_option('output_filename')  # or get_attribute(), but this is clearer
            with retrieved.open(filename_stdout, 'r') as fil:
                out_file = fil.readlines()
        except OSError:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        job_done = False
        for i in range(len(out_file)):
            line = out_file[-i]
            if 'JOB DONE' in line:
                job_done = True
                break
        if not job_done:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE)

        # Parse basic info and warnings, and output them as output_parmeters
        parsed_data, logs = parse_output_base(out_file, 'PROJWFC')
        self.emit_logs(logs)
        self.out('output_parameters', Dict(dict=parsed_data))

        # check and read pdos_tot file
        out_filenames = retrieved.list_object_names()
        try:
            pdostot_filename = fnmatch.filter(out_filenames, '*pdos_tot*')[0]
            with retrieved.open(pdostot_filename, 'r') as pdostot_file:
                # Columns: Energy(eV), Ldos, Pdos
                pdostot_array = np.atleast_2d(np.genfromtxt(pdostot_file))
                energy = pdostot_array[:, 0]
                dos = pdostot_array[:, 1]
        except (OSError, KeyError):
            return self.exit(self.exit_codes.ERROR_READING_PDOSTOT_FILE)

        # check and read all of the individual pdos_atm files
        pdos_atm_filenames = fnmatch.filter(out_filenames, '*pdos_atm*')
        pdos_atm_array_dict = {}
        for name in pdos_atm_filenames:
            with retrieved.open(name, 'r') as pdosatm_file:
                pdos_atm_array_dict[name] = np.atleast_2d(np.genfromtxt(pdosatm_file))

        # finding the bands and projections
        # we create a dictionary the progressively accumulates more info
        out_info_dict = {}
        out_info_dict['out_file'] = out_file
        out_info_dict['energy'] = energy
        out_info_dict['pdos_atm_array_dict'] = pdos_atm_array_dict
        try:
            new_nodes_list = self._parse_bands_and_projections(out_info_dict)
        except QEOutputParsingError as err:
            self.logger.error(f'Error parsing bands and projections: {err}')
            return self.exit(self.exit_codes.ERROR_PARSING_PROJECTIONS)
        for linkname, node in new_nodes_list:
            self.out(linkname, node)

        Dos_out = XyData()
        Dos_out.set_x(energy, 'Energy', 'eV')
        Dos_out.set_y(dos, 'Dos', 'states/eV')
        self.out('Dos', Dos_out)

    def _parse_bands_and_projections(self, out_info_dict):
        """Function that parses the standard output into bands and projection data.

        :param out_info_dict: used to pass technical internal variables
                              to helper functions in compact form
        :return: append_nodes_list a list containing BandsData and
                 ProjectionData parsed from standard_out
        """
        out_file = out_info_dict['out_file']  # Note: we expect a list of lines
        out_info_dict['k_lines'] = []
        out_info_dict['e_lines'] = []
        out_info_dict['psi_lines'] = []
        out_info_dict['wfc_lines'] = []
        append_nodes_list = []

        for i, line in enumerate(out_file):
            if 'k =' in line:
                out_info_dict['k_lines'].append(i)
            if '==== e(' in line or line.strip().startswith('e ='):
                # The energy format in output was changed in QE6.3
                # this check supports old and new format
                out_info_dict['e_lines'].append(i)
            if '|psi|^2' in line:
                out_info_dict['psi_lines'].append(i)
            if 'state #' in line:
                out_info_dict['wfc_lines'].append(i)

        # Basic check
        if len(out_info_dict['e_lines']) != len(out_info_dict['psi_lines']):
            raise QEOutputParsingError('e-lines and psi-lines are in different number')
        if len(out_info_dict['psi_lines']) % len(out_info_dict['k_lines']) != 0:
            raise QEOutputParsingError('Band Energy Points is not a multiple of kpoints')
        # calculates the number of bands
        out_info_dict['num_bands'] = len(out_info_dict['psi_lines']) // len(out_info_dict['k_lines'])

        # Uses the parent input parameters, and checks if the parent used
        # spin calculations. Try to replace with a query, if possible.
        try:
            parent_calc = (
                self.node.inputs.parent_folder.get_incoming(node_class=CalcJobNode,
                                                            link_type=LinkType.CREATE).one().node
            )
        except ValueError as e:
            raise QEOutputParsingError(f'Could not get parent calculation of input folder: {e}')
        out_info_dict['parent_calc'] = parent_calc
        try:
            parent_param = parent_calc.get_outgoing(link_label_filter='output_parameters').one().node
        except ValueError:
            raise QEOutputParsingError('The parent had no output_parameters! Cannot parse from this!')
        try:
            structure = parent_calc.get_incoming(link_label_filter='structure').one().node
        except ValueError:
            raise QEOutputParsingError('The parent had no input structure! Cannot parse from this!')
        try:
            nspin = parent_param.get_dict()['number_of_spin_components']
            if nspin != 1:
                spin = True
            else:
                spin = False
            out_info_dict['spinorbit'] = parent_param.get_dict().get('spin_orbit_calculation', False)
            out_info_dict['collinear'] = not parent_param.get_dict().get('non_colinear_calculation', False)
            if not out_info_dict['collinear']:
                # Sanity check
                if nspin != 4:
                    raise QEOutputParsingError('The calculation is non-collinear, but nspin is not set to 4!')
                spin = False
        except KeyError:
            spin = False
            out_info_dict['spinorbit'] = False
            out_info_dict['collinear'] = True
        out_info_dict['spin'] = spin

        # changes k-numbers to match spin
        # because if spin is on, k points double for up and down
        out_info_dict['k_states'] = len(out_info_dict['k_lines'])
        if spin:
            if out_info_dict['k_states'] % 2 != 0:
                raise QEOutputParsingError('Internal formatting error regarding spin')
            out_info_dict['k_states'] = out_info_dict['k_states'] // 2

        #   adds in the k-vector for each kpoint
        k_vect = [out_file[out_info_dict['k_lines'][i]].split()[2:] for i in range(out_info_dict['k_states'])]
        out_info_dict['k_vect'] = np.array(k_vect)
        out_info_dict['structure'] = structure
        out_info_dict['orbitals'] = find_orbitals_from_statelines(out_info_dict)

        if spin:
            # I had to guess what the ordering of the spin is, because
            # the projwfc.x documentation doesn't say, but looking at the
            # source code I found:
            #
            # DO is=1,nspin
            #   IF (nspin==2) THEN
            #       IF (is==1) filename=trim(filproj)//'.up'
            #       IF (is==2) filename=trim(filproj)//'.down'
            #
            # Which would say that it is reasonable to assume that the
            # spin up states are written first, then spin down
            #
            out_info_dict['spin_down'] = False
            bands_data1, projection_data1 = spin_dependent_subparser(out_info_dict)
            append_nodes_list += [('projections_up', projection_data1), ('bands_up', bands_data1)]
            out_info_dict['spin_down'] = True
            bands_data2, projection_data2 = spin_dependent_subparser(out_info_dict)
            append_nodes_list += [('projections_down', projection_data2), ('bands_down', bands_data2)]
        else:
            out_info_dict['spin_down'] = False
            bands_data, projection_data = spin_dependent_subparser(out_info_dict)
            append_nodes_list += [('projections', projection_data), ('bands', bands_data)]

        return append_nodes_list
