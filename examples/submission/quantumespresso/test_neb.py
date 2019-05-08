#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
###########################################################################
# Copyright (c), The AiiDA team. All rights reserved.                     #
# This file is part of the AiiDA code.                                    #
#                                                                         #
# The code is hosted on GitHub at https://github.com/aiidateam/aiida_core #
# For further information on the license, see the LICENSE.txt file        #
# For further information please visit http://www.aiida.net               #
###########################################################################
from __future__ import absolute_import
from __future__ import print_function
import sys
import os

from aiida.common.example_helpers import test_and_get_code

from aiida.common.exceptions import NotExistent

################################################################

UpfData = DataFactory('upf')
Dict = DataFactory('dict')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')
try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print(("The first parameter can only be either "
                          "--send or --dont-send"), file=sys.stderr)
    sys.exit(1)

try:
    codename = sys.argv[2]
except IndexError:
    codename = None

# If True, load the pseudos from the family specified below
# Otherwise, use static files provided
auto_pseudos = True

queue = None
# queue = "Q_aries_free"
settings = None
#####

code = test_and_get_code(codename, expected_code_type='quantumespresso.neb')

cell = [[4.0, 0., 0., ],
        [0., 4.0, 0., ],
        [0., 0.,  4.05, ],]

# Displacements
d1 = 1.940525216
d2 = 0.183506832
d3 = 2.216318307

# Initial structure
s1 = StructureData(cell=cell)
s1.append_atom(position=(0., 0., 0.), symbols=['Ba'])
s1.append_atom(position=(2., 2., d1), symbols=['Ti'])
s1.append_atom(position=(2., 2., d2), symbols=['O'])
s1.append_atom(position=(2., 0., d3), symbols=['O'])
s1.append_atom(position=(0., 2., d3), symbols=['O'])

# Final structure
s2 = StructureData(cell=cell)
s2.append_atom(position=(0., 0., 0.), symbols=['Ba'])
s2.append_atom(position=( 2., 2., 4.05 - d1), symbols=['Ti'])
s2.append_atom(position=( 2., 2., -d2), symbols=['O'])
s2.append_atom(position=( 2., 0., 4.05 - d3), symbols=['O'])
s2.append_atom(position=( 0., 2., 4.05 - d3), symbols=['O'])

fixed_coords = [] 
for site in s1.sites:
    if site.kind_name == 'Ba':
        fixed_coords.append([True,True,True])
    else:
        fixed_coords.append([True,True,False])
settings = Dict(dict={
        'fixed_coords': fixed_coords
        })

elements = list(s1.get_symbols_set())

valid_pseudo_groups = UpfData.get_upf_groups(filter_elements=elements)


if auto_pseudos:
    valid_pseudo_groups = UpfData.get_upf_groups(filter_elements=elements)

    try:
        pseudo_family = sys.argv[3]
    except IndexError:
        print("Error, auto_pseudos set to True. You therefore need to pass as second parameter", file=sys.stderr)
        print("the pseudo family name.", file=sys.stderr)
        print("Valid UPF families are:", file=sys.stderr)
        print("\n".join("* {}".format(i.name) for i in valid_pseudo_groups), file=sys.stderr)
        sys.exit(1)

    try:
        UpfData.get_upf_group(pseudo_family)
    except NotExistent:
        print("auto_pseudos is set to True and pseudo_family='{}',".format(pseudo_family), file=sys.stderr)
        print("but no group with such a name found in the DB.", file=sys.stderr)
        print("Valid UPF groups are:", file=sys.stderr)
        print(",".join(i.name for i in valid_pseudo_groups), file=sys.stderr)
        sys.exit(1)
    
pw_parameters = Dict(dict={
    'CONTROL': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
    },
    'SYSTEM': {
        'ecutwfc': 30.,
        'ecutrho': 240.,
    },
    'ELECTRONS': {
        'conv_thr': 1.e-8,
        'mixing_beta': 0.3,

    }})

neb_parameters = Dict(dict={
        'PATH': {
            'restart_mode': 'from_scratch',
            'string_method': 'neb',
            'nstep_path': 20, 
            'ds': 2.0,
            'opt_scheme': 'broyden',
            'num_of_images': 7,
            'k_max': 0.3,
            'k_min': 0.2,
            'path_thr': 0.2,
        }})

# If you want to set a manual climbing image:
#neb_parameters.dict.PATH['ci_scheme']=  'manual'
# and specify the climbing image(s) index in settings
#try:
#    settings.update_dict({'climbing_images': [4]})
#except NameError:
#    settings = Dict(dict={'climbing_images': [4]})

kpoints = KpointsData()

# method mesh
kpoints_mesh = 2
kpoints.set_kpoints_mesh([kpoints_mesh, kpoints_mesh, kpoints_mesh])


calc = code.new_calc()
calc.label = "Test QE neb.x"
calc.description = "Test calculation with the Quantum ESPRESSO neb.x code"
calc.set_option('max_wallclock_seconds', 30 * 60)  # 30 min
# Valid only for Slurm and PBS (using default values for the
# number_cpus_per_machine), change for SGE-like schedulers 
calc.set_option('resources', {"num_machines": 1})
## Otherwise, to specify a given # of cpus per machine, uncomment the following:
# calc.set_option('resources', {"num_machines": 1, "num_mpiprocs_per_machine": 8})

#calc.set_option('custom_scheduler_commands', ("#SBATCH --account=ch3")

if queue is not None:
    calc.set_option('queue_name', queue)

calc.use_first_structure(s1)
calc.use_last_structure(s2)
calc.use_pw_parameters(pw_parameters)
calc.use_neb_parameters(neb_parameters)

try:
    calc.use_pseudos_from_family(pseudo_family)
    print("Pseudos successfully loaded from family {}".format(pseudo_family))
except NotExistent:
    print ("Pseudo or pseudo family not found. You may want to load the "
           "pseudo family, or set auto_pseudos to False.")
    raise

calc.use_kpoints(kpoints)

if settings is not None:
    calc.use_settings(settings)
#from aiida.orm.nodes.data.remote import RemoteData
#calc.set_outdir(remotedata)

if submit_test:
    subfolder, script_filename = calc.submit_test()
    print("Test_submit for calculation (uuid='{}')".format(
        calc.uuid))
    print("Submit file in {}".format(os.path.join(
        os.path.relpath(subfolder.abspath),
        script_filename
    )))
else:
    calc.store_all()
    print("created calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid, calc.dbnode.pk))
    calc.submit()
    print("submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid, calc.dbnode.pk))

