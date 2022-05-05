#!/usr/bin/env python
# -*- coding: utf-8 -*-
###########################################################################
# Copyright (c), The AiiDA team. All rights reserved.                     #
# This file is part of the AiiDA code.                                    #
#                                                                         #
# The code is hosted on GitHub at https://github.com/aiidateam/aiida_core #
# For further information on the license, see the LICENSE.txt file        #
# For further information please visit http://www.aiida.net               #
###########################################################################
# start-marker for docs

from aiida import engine, load_profile, orm, plugins

load_profile()

builder = orm.load_code('pw-6.3@TheHive').get_builder()

# BaTiO3 cubic structure
alat = 4.  # angstrom
cell = [[alat, 0., 0.], [0., alat, 0.], [0., 0., alat]]
s = plugins.DataFactory('core.structure')(cell=cell)
s.append_atom(position=(0., 0., 0.), symbols='Ba')
s.append_atom(position=(alat / 2., alat / 2., alat / 2.), symbols='Ti')
s.append_atom(position=(alat / 2., alat / 2., 0.), symbols='O')
s.append_atom(position=(alat / 2., 0., alat / 2.), symbols='O')
s.append_atom(position=(0., alat / 2., alat / 2.), symbols='O')
builder.structure = s
builder.pseudos = orm.load_group('SSSP/1.1/PBE/efficiency').get_pseudos(structure=s)

builder.parameters = plugins.DataFactory('core.dict')(
    dict={
        'CONTROL': {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'wf_collect': True,
        },
        'SYSTEM': {
            'ecutwfc': 30.,
            'ecutrho': 240.,
        },
        'ELECTRONS': {
            'conv_thr': 1.e-6,
        }
    }
)

kpoints = plugins.DataFactory('core.array.kpoints')()
kpoints.set_kpoints_mesh([4, 4, 4])
builder.kpoints = kpoints

builder.metadata.label = 'BaTiO3 test run'
builder.metadata.options.resources = {'num_machines': 1}
builder.metadata.options.max_wallclock_seconds = 1800

calc = engine.submit(builder)
print(f'created calculation with PK={calc.pk}')
