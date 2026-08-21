# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the d3hess.x code of Quantum ESPRESSO."""

from aiida import orm

from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class D3hessCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the d3hess.x code of Quantum ESPRESSO.

    d3hess.x code of the Quantum ESPRESSO distribution computes the Hessian matrix of the DFT-D3 dispersion
    after an SCF calculation.
    """

    _default_namelists = ['INPUT']
    _blocked_keywords = [
        ('INPUT', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),  # pylint: disable=protected-access
        ('INPUT', 'prefix', NamelistsCalculation._PREFIX),  # pylint: disable=protected-access
    ]

    _internal_retrieve_list = []

    _DEFAULT_INPUT_FILE = 'd3hess.in'
    _DEFAULT_OUTPUT_FILE = 'd3hess.out'

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True)
        # With symlinking, the file `aiida.hess` produced in ./out/ will also be created in the parent folder.
        # This allows to reuse the remote folder of the SCF calculation (instead of the newly created one for d3hess)
        # for post-processings which require the `aiida.hess` file (e.g. ph.x with dft-d3 correction)
        # This prevents the eventual unintended killing of the logic of post-processing workflows which can depend on
        # the parenthood of the remote folder to an actual SCF calculation
        spec.input('settings', valid_type=orm.Dict, default=lambda: orm.Dict(dict={'PARENT_FOLDER_SYMLINK': True}))
        spec.input('metadata.options.withmpi', valid_type=bool, default=False)  # Usually takes less than a sec
        # yapf: enable
