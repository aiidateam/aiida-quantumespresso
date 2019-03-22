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
import sys, os, numpy
from aiida.common.example_helpers import test_and_get_code

Dict = DataFactory('dict')

try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print >> sys.stderr, ("The first parameter can only be either "
                          "--send or --dont-send")
    sys.exit(1)

try:
    parent_id = sys.argv[2]
    nnkp_file_id = sys.argv[3]
    codename = sys.argv[4]
except IndexError:
    print >> sys.stderr, ("Must provide as further parameters:\n- the parent ID"
        "\n- the nnkp SingleFile\n- "
        "a pw2wannier90.x codename")
    sys.exit(1)

num_machines = 1 # node numbers

#####
try:
    parent_id = int(parent_id)
except ValueError:
    print >> sys.stderr, 'Parent_id not an integer: {}'.format(parent_id)
    sys.exit(1)

try:
    nnkp_file_id = int(nnkp_file_id)
except ValueError:
    print >> sys.stderr, 'nnkp_file_id not an integer: {}'.format(nnkp_file_id)
    sys.exit(1)

nnkp_file = load_node(nnkp_file_id)
if not isinstance(nnkp_file, DataFactory('singlefile')):
    print >> sys.stderr, 'The provided nnkp_file is not an SinglefileData: {}'.format(nnkp_file_id)
    sys.exit(1)    

code = test_and_get_code(codename, expected_code_type='quantumespresso.pw2wannier90')

computer = code.get_remote_computer()

parameters = Dict(dict={
        'INPUTPP': {
            'write_unk': True,
            'write_amn': True,
            'write_mmn': True,
        },
     })

parentcalc = load_node(parent_id)

# calc = code.new_calc(computer=computer)
calc = code.new_calc()
calc.label = "Test QE pw2wannier90.x"
calc.description = "Test calculation with the Quantum ESPRESSO pw2wannier90.x code"
calc.set_option('max_wallclock_seconds', 60*30) # 30 min
calc.set_option('resources', {"num_machines":num_machines})

calc.use_parameters(parameters)
calc.use_parent_calculation(parentcalc)
calc.use_nnkp_file(nnkp_file)

if submit_test:
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(
        calc.uuid)
    print "Submit file in {}".format(os.path.join(
        os.path.relpath(subfolder.abspath),
        script_filename
        ))
else:
    calc.store_all()
    print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.pk)
    calc.submit()
    print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.pk)
