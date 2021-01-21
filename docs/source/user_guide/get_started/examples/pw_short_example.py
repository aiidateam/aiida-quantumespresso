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
from aiida import load_profile

from aiida.orm import Code
from aiida.plugins import DataFactory
from aiida.engine import submit
from aiida.orm.nodes.data.upf import get_pseudos_from_structure

load_profile()

StructureData = DataFactory('structure')
Dict = DataFactory('dict')
KpointsData = DataFactory('array.kpoints')

###############################
# Set your values here
codename = 'pw-6.3@TheHive'
pseudo_family = 'SSSP_efficiency'
###############################

code = Code.get_from_string(codename)
builder = code.get_builder()

# BaTiO3 cubic structure
alat = 4.  # angstrom
cell = [
    [
        alat,
        0.,
        0.,
    ],
    [
        0.,
        alat,
        0.,
    ],
    [
        0.,
        0.,
        alat,
    ],
]
s = StructureData(cell=cell)
s.append_atom(position=(0., 0., 0.), symbols='Ba')
s.append_atom(position=(alat / 2., alat / 2., alat / 2.), symbols='Ti')
s.append_atom(position=(alat / 2., alat / 2., 0.), symbols='O')
s.append_atom(position=(alat / 2., 0., alat / 2.), symbols='O')
s.append_atom(position=(0., alat / 2., alat / 2.), symbols='O')

parameters = Dict(
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

kpoints = KpointsData()
kpoints.set_kpoints_mesh([4, 4, 4])

builder.pseudos = get_pseudos_from_structure(s, pseudo_family)

builder.metadata.options.resources = {'num_machines': 1}
builder.metadata.options.max_wallclock_seconds = 1800

builder.metadata.label = 'My generic title'
builder.metadata.description = 'My generic description'

builder.structure = s
builder.parameters = parameters
builder.kpoints = kpoints

calc = submit(builder)

print(f'created calculation with PK={calc.pk}')
