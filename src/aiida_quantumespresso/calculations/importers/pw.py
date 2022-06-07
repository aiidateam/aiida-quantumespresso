# -*- coding: utf-8 -*-
""":class:`aiida.engine.processes.calcjobs.importer.CalcJobImporter` implementation for ``PwCalculation``."""
from __future__ import annotations

import pathlib
import tempfile

from aiida.engine import CalcJobImporter
from aiida.orm import AbstractCode, Node, RemoteData


class PwCalculationImporter(CalcJobImporter):
    """:class:`aiida.engine.processes.calcjobs.importer.CalcJobImporter` implementation for ``PwCalculation``.

    This class allows to import a completed ``pw.x`` calculation that was executed without AiiDA, into an AiiDA profile.
    """

    @staticmethod
    def parse_remote_data(  # pylint: disable=arguments-differ
        remote_data: RemoteData,
        input_file_name: str,
        code: AbstractCode | None = None,
        metadata: dict | None = None,
        pseudo_folder_path: str | None = None,
    ) -> dict[str, Node | dict]:
        """Parse the input nodes from the files in the provided ``RemoteData``.

        :param remote_data: the remote data node containing the raw input files.
        :param input_file_name: the filename of the main Quantum ESPRESSO input file.
        :param code: Optional ``AbstractCode`` node to attach to the imported node, as if it would have been run with
            the ``code`` input in a real run.
        :param metadata: Optional ``metadata`` inputs to set on the imported node, as if it would have been run with
            the ``metadata`` input in a real run.
        :returns: a dictionary with the parsed inputs nodes that match the input spec of the associated ``CalcJob``.
        """
        from aiida_quantumespresso.tools.pwinputparser import create_builder_from_file

        with tempfile.TemporaryDirectory() as tmppath:
            dirpath = pathlib.Path(tmppath) / 'folder'
            with remote_data.get_authinfo().get_transport() as transport:
                transport.copytree(remote_data.get_remote_path(), dirpath)

            builder = create_builder_from_file(
                str(dirpath), input_file_name, code, metadata or {}, pseudo_folder_path=pseudo_folder_path
            )

        inputs = dict(builder)
        inputs['remote_folder'] = remote_data

        return inputs
