# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the bands.x code of Quantum ESPRESSO."""

from aiida import orm

from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class BandsCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the bands.x code of Quantum ESPRESSO.

    bands.x code of the Quantum ESPRESSO distribution, re-orders bands, and computes band-related properties.

    It computes for instance the expectation value of the momentum operator:
    <Psi(n,k) | i * m * [H, x] | Psi(m,k)>. For more information, refer to http://www.quantum-espresso.org/
    """

    _MOMENTUM_OPERATOR_NAME = 'momentum_operator.dat'
    _BANDS_NAME = 'bands.dat'

    _default_namelists = ['BANDS']
    _blocked_keywords = [
        ('BANDS', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),  # pylint: disable=protected-access
        ('BANDS', 'prefix', NamelistsCalculation._PREFIX),  # pylint: disable=protected-access
        ('BANDS', 'filband', _BANDS_NAME),
        ('BANDS', 'filp', _MOMENTUM_OPERATOR_NAME),  # Momentum operator
    ]

    _internal_retrieve_list = []
    _default_parser = 'quantumespresso.bands'

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True)
        spec.output('output_parameters', valid_type=orm.Dict)
        # yapf: enable
