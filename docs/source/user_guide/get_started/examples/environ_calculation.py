# -*- coding: utf-8 -*-

from aiida import orm
from aiida.engine import run, submit
import numpy as np
from ase import Atoms
from ase import units
"""
Calculation that recreates the environ calculation of
CO on Pt(111) with a unit charge

More details about the solvation models and potential corrections
can be found here:

O. Andreussi, I. Dabo and N. Marzari, J. Chem. Phys. 136,064102 (2012).
I. Dabo et al. Phys.Rev. B 77, 115139 (2008)

This system is taken directly from Example04 in the
examples folder from Environ
"""


def calculator():
    """
    This is the SCF calculator, which is identical
    and unchanged from what one would do in the case
    of a standard SCF calculation - nothing environ
    related. Remember to add in the tot_charge here
    """
    param_dict = {
        'CONTROL': {
            'calculation': 'scf',
            'tprnfor': True,
        },
        'SYSTEM': {
            'ecutwfc': 35,
            'ecutrho': 280,
            'occupations': 'smearing',
            'smearing': 'cold',
            'degauss': 0.03,
            'tot_charge': 1,
        },
        'ELECTRONS': {
            'conv_thr': 5e-9,
            'electron_maxstep': 200,
            'mixing_beta': 0.2,
            'mixing_ndim': 15,
            'diagonalization': 'david',
        },
    }

    return param_dict


def environ_calculator():
    """
    Environ parameters
    """
    param_dict = {
        'ENVIRON': {
            'verbose': 0,
            'environ_thr': 1.0,
            'environ_type': 'input',
            'env_static_permittivity': 1,
            'env_surface_tension': 0.0,
            'env_pressure': 0.0,
            'env_electrostatic': True,
        },
        'BOUNDARY': {
            'solvent_mode': 'full',
        },
        'ELECTROSTATIC': {
            'pbc_correction': 'parabolic',
            'pbc_dim': 2,
            'pbc_axis': 3,
            'tol': 5e-13,
        }
    }
    return param_dict


def runner(code, structure):

    ## Base calculation for testing
    PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
    builder = PwBaseWorkChain.get_builder()
    builder.pw.structure = structure

    builder.metadata.label = 'Single point environ'
    builder.metadata.description = 'Testing calculation with environ for Pt(111) + CO'

    KpointsData = DataFactory('array.kpoints')
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([1, 1, 1])  # Not physical just a test
    builder.kpoints = kpoints

    family = load_group('SSSP/1.1/PBE/efficiency')
    builder.pw.pseudos = family.get_pseudos(structure=structure)

    calculation = calculator()
    environ_calculation = environ_calculator()
    ## these are the main cube files that could potentially be parsed
    ## if the verbosity is set to 2 or higher
    # environ_calculation['additional_retrieve_list'] = ['epsilon.cube', \
    #     'vreference.cube', 'velectrostatic.cube', 'vsoftcavity.cube', \
    #     'electrons.cube', 'charges.cube', 'smeared_ions.cube']

    builder.pw.parameters = orm.Dict(dict=calculation)
    builder.pw.settings = orm.Dict(dict=environ_calculation)
    builder.pw.metadata.options.resources = {'num_machines': 1}
    builder.pw.metadata.options.max_wallclock_seconds = 5 * 60
    builder.pw.code = code

    calculation = submit(builder)


if __name__ == '__main__':
    code = load_code('pw_6-7@juwels')

    StructureData = DataFactory('structure')
    ## these are the original coordinates for the Pt-CO system
    positions = [
        [5.335084148, 4.646723426, 12.901029877],
        [5.335009643, 4.619623254, 15.079854269],
        [8.061327071, 0.098057998, 8.992142901],
        [2.608989366, 0.098058283, 8.992140585],
        [0.000036609, 4.720846294, 8.968756935],
        [5.335159557, 4.721612729, 9.380196435],
        [0.000041121, 7.802951963, 4.604626508],
        [5.335161233, 7.697749113, 4.753489408],
        [2.697860636, 3.152173889, 4.688412329],
        [7.972463687, 3.152174491, 4.688415209],
    ]

    ## setting up the system with ASE
    ## notice the units that are being used
    atoms = Atoms('COPt8')
    atoms.set_positions(np.array(positions) * units.Bohr)
    atoms.set_pbc([True, True, True])
    a = 10.6881 * units.Bohr
    b = 0.866025 * a * units.Bohr
    c = 3.95422 * a * units.Bohr
    atoms.set_cell([a, b, c])
    structure = StructureData(ase=atoms)

    runner(code=code, structure=structure)
