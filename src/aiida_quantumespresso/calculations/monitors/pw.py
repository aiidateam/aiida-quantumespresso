# -*- coding: utf-8 -*-
"""Monitors for :class:`aiida_quantumespresso.calculations.pw.PwCalculation` jobs."""
from __future__ import annotations

import re
import tempfile

from aiida.orm import CalcJobNode
from aiida.transports import Transport


def monitor_scf_iterations(
    node: CalcJobNode, transport: Transport, max_number_iterations: int | None = None
) -> str | None:
    """Monitor the number of SCF iterations and interrupt when it exceeds ``max_number_iterations``.

    Note that this monitor is just for demonstration purposes. If you want to limit the number of SCF convergence
    iterations, it is better to use the ``electron_maxstep`` keyword in the ``ELECTRONS`` card of the ``parameters``
    input.

    :param max_number_iterations: If the number of SCF convergence steps exceeds this number, kill the job.
    """
    with tempfile.NamedTemporaryFile('w+') as handle:
        transport.getfile(node.get_option('output_filename'), handle.name)
        handle.seek(0)
        output = handle.read()

    if max_number_iterations is not None and re.search(fr'iteration #\s*{max_number_iterations}', output):
        return f'Maximum number of iterations {max_number_iterations} in SCF convergence cycle reached.'
