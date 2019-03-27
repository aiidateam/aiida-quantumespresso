# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.common.extendeddicts import AttributeDict

pw = AttributeDict({
    'conv_thr': 1e-6,
    'degauss': 0.,
    'diagonalization': 'david',
    'electron_maxstep': 100,
    'mixing_beta': 0.7,
    'mixing_mode': 'plain',
    'mixing_ndim': 8,
    'noncolin': False,
    'nspin': 1,
    'occupations': None,
    'press': 0.,
    'press_conv_thr': 0.5,
    'smearing': '',
    'startmag': 0.,
    'wf_collect': False,
})