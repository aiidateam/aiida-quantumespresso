#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import time

from aiida.orm import DataFactory

# If set to True, will ask AiiDA to run in serial mode (i.e., AiiDA will not
# invoke the mpirun command in the submission script)
run_in_serial_mode = False
codename = 'qe-pw@torquessh'
# If it takes > 5 min, I decide I failed (e.g., no daemon is running)
timeout_secs = 5*60 
queue = None

expected_energy = -3700.91106342615

################################################################

UpfData = DataFactory('upf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

code = Code.get_from_string(codename)

alat = 4.  # angstrom
cell = [[alat, 0., 0., ],
        [0., alat, 0., ],
        [0., 0., alat, ],
]

# BaTiO3 cubic structure
s = StructureData(cell=cell)
s.append_atom(position=(0., 0., 0.), symbols=['Ba'])
s.append_atom(position=(alat / 2., alat / 2., alat / 2.), symbols=['Ti'])
s.append_atom(position=(alat / 2., alat / 2., 0.), symbols=['O'])
s.append_atom(position=(alat / 2., 0., alat / 2.), symbols=['O'])
s.append_atom(position=(0., alat / 2., alat / 2.), symbols=['O'])

elements = list(s.get_symbols_set())

parameters = ParameterData(dict={
    'CONTROL': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'wf_collect': True,
        'tstress': True,
        'tprnfor': True,
    },
    'SYSTEM': {
        'ecutwfc': 40.,
        'ecutrho': 320.,
    },
    'ELECTRONS': {
        'conv_thr': 1.e-10,
    }})

kpoints = KpointsData()
kpoints_mesh = 2
kpoints.set_kpoints_mesh([kpoints_mesh, kpoints_mesh, kpoints_mesh])

# to retrieve the bands
# (the object settings is optional)
settings_dict = {}
settings = ParameterData(dict=settings_dict)

calc = code.new_calc()
calc.label = "Test QE pw.x"
calc.description = "Test calculation with the Quantum ESPRESSO pw.x code"
calc.set_max_wallclock_seconds(30 * 60)  # 30 min
# Valid only for Slurm and PBS (using default values for the
# number_cpus_per_machine), change for SGE-like schedulers
calc.set_resources({"num_machines": 1})
if run_in_serial_mode:
    calc.set_withmpi(False)

if queue is not None:
    calc.set_queue_name(queue)

calc.use_structure(s)
calc.use_parameters(parameters)

raw_pseudos = [
    ("Ba.pbesol-spn-rrkjus_psl.0.2.3-tot-pslib030.UPF", 'Ba', 'pbesol'),
    ("Ti.pbesol-spn-rrkjus_psl.0.2.3-tot-pslib030.UPF", 'Ti', 'pbesol'),
    ("O.pbesol-n-rrkjus_psl.0.1-tested-pslib030.UPF", 'O', 'pbesol')]

pseudos_to_use = {}
for fname, elem, pot_type in raw_pseudos:
    absname = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                            "data", fname))
    pseudo, created = UpfData.get_or_create(
        absname, use_first=True)
    if created:
        print "Created the pseudo for {}".format(elem)
    else:
        print "Using the pseudo for {} from DB: {}".format(elem, pseudo.pk)
    pseudos_to_use[elem] = pseudo

for k, v in pseudos_to_use.iteritems():
    calc.use_pseudo(v, kind=k)

calc.use_kpoints(kpoints)

if settings is not None:
    calc.use_settings(settings)

calc.store_all()
print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
    calc.uuid, calc.dbnode.pk)
calc.submit()
print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
    calc.uuid, calc.dbnode.pk)


print "Wating for end of execution..."
start_time = time.time()
exited_with_timeout = True
while time.time() - start_time < timeout_secs:
    time.sleep(15) # Wait a few seconds
    

    # print some debug info, both for debugging reasons and to avoid
    # that the test machine is shut down because there is no output

    print "#"*78
    print "####### TIME ELAPSED: {} s".format(time.time() - start_time)
    print "#"*78
    print "Output of 'verdi calculation list':"
    try:
        print subprocess.check_output(
            ["verdi", "calculation", "list"], 
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print "Note: the command failed, message: {}".format(e.message)

    if calc.is_terminated:
        print "Calculation terminated its execution"
        exited_with_timeout = False
        break

if exited_with_timeout:
    print "Timeout!! Calculation did not complete after {} seconds".format(
        timeout_secs)
    sys.exit(2)
else:
    if abs(calc.res.energy - expected_energy) < 1.e-3:
        print "OK, energy has the expected value"
        sys.exit(0)
    else:
        print "ERROR!"
        print "Expected energy value: {}".format(expected_energy)
        print "Actual energy value: {}".format(calc.res.energy)
        sys.exit(3)
        
        
