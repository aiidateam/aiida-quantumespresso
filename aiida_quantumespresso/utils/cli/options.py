# -*- coding: utf-8 -*-
from aiida.utils.cli.options import overridable_option

automatic_parallelization = overridable_option(
    '-a', '--automatic-parallelization', is_flag=True, default=False, show_default=True,
    help='enable the automatic parallelization option of the workchain'
)

clean_workdir = overridable_option(
    '-x', '--clean-workdir', is_flag=True, default=False, show_default=True,
    help='clean the remote folder of all the launched calculations after completion of the workchain'
)