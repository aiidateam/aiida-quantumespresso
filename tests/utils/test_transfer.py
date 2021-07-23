# -*- coding: utf-8 -*-
"""Unit tests for the :py:mod:`~aiida_quantumespresso.utils.transfer` module."""
import pytest
from packaging import version
import aiida

from aiida.orm import FolderData, RemoteData
from aiida.engine import ProcessBuilder
from aiida_quantumespresso.utils.transfer import get_transfer_builder

pytestmark = pytest.mark.skipif(
    version.parse(aiida.get_version()) < version.parse('1.6'),
    reason='Transfer was released in AiiDA core v1.6',
)


@pytest.fixture(name='make_data_source')
def fixture_make_data_source(make_tempdir):
    """Provide a function to generate data sources specific for this test."""

    def _make_data_source(computer=None, populate_level=0):
        import os
        import io

        if computer is None:
            node = FolderData()

            filepaths = []
            if populate_level > 0:
                filepaths.append('data-file-schema.xml')
                filepaths.append('charge-density.dat')
            if populate_level > 1:
                filepaths.append('paw.txt')

            for filepath in filepaths:
                node.put_object_from_filelike(io.StringIO(), path=filepath)

        else:
            remote_path = make_tempdir(parent_dir=computer.get_workdir())
            node = RemoteData(computer=computer, remote_path=remote_path)

            filepaths = []
            if populate_level > 0:
                filepaths.append('out/aiida.save/data-file-schema.xml')
                filepaths.append('out/aiida.save/charge-density.dat')
            if populate_level > 1:
                filepaths.append('out/aiida.save/paw.txt')

            for filepath in filepaths:
                fullpath = os.path.join(remote_path, filepath)
                os.makedirs(os.path.dirname(fullpath), exist_ok=True)
                open(fullpath, 'w').close()

        return node.store()

    return _make_data_source


@pytest.mark.parametrize('track', [True, False])
@pytest.mark.parametrize('populate_level', [1, 2])
@pytest.mark.parametrize('source_is_remote, provide_computer', [(True, True), (True, False), (False, True)])
def test_transfer_builder_results(
    makeget_computer, make_data_source, source_is_remote, provide_computer, populate_level, track
):
    """Test all viable input sets for the get_transfer_builder utility function."""

    if provide_computer:
        builder_computer = makeget_computer('computer1')
        check_computer = builder_computer
    else:
        builder_computer = None

    if source_is_remote:
        folder_computer = makeget_computer('computer2')
        data_source = make_data_source(computer=folder_computer, populate_level=populate_level)
        retrieve_expected = True
        expected_listname = 'symlink_files'
        check_computer = folder_computer
    else:
        data_source = make_data_source(computer=None, populate_level=populate_level)
        retrieve_expected = False
        expected_listname = 'local_files'

    if source_is_remote and provide_computer:
        with pytest.warns(UserWarning) as warnings:
            builder = get_transfer_builder(data_source, computer=builder_computer, track=track)
        assert len(warnings) == 1
        assert 'ignore' in str(warnings[0].message)
        assert f'{builder_computer}' in str(warnings[0].message)
        assert f'{folder_computer}' in str(warnings[0].message)
        # Makes sure the information is there, but there is no generic way of checking the correctedness
        # of the warning (i.e. which computer is the one selected) without checking for the exact text
    else:
        builder = get_transfer_builder(data_source, computer=builder_computer, track=track)

    assert isinstance(builder, ProcessBuilder)
    assert builder.metadata['computer'].pk == check_computer.pk
    assert builder.source_nodes['source_node'] == data_source
    assert builder.instructions.get_dict()['retrieve_files'] == retrieve_expected
    assert expected_listname in builder.instructions.get_dict()

    if source_is_remote:
        # paw.txt are always loaded when remote copying
        expected_setlist = {
            ('source_node', 'out/aiida.save/data-file-schema.xml', 'data-file-schema.xml'),
            ('source_node', 'out/aiida.save/charge-density.dat', 'charge-density.dat'),
            ('source_node', 'out/aiida.save/paw.txt', 'paw.txt'),
        }

    else:
        expected_setlist = {
            ('source_node', 'data-file-schema.xml', 'out/aiida.save/data-file-schema.xml'),
            ('source_node', 'charge-density.dat', 'out/aiida.save/charge-density.dat'),
        }
        if populate_level == 2:
            expected_setlist.add(('source_node', 'paw.txt', 'out/aiida.save/paw.txt'))

    # Note: I need to do the following transformation manually because sometimes what I get from:
    #
    # > builder.instructions.get_dict()[expected_listname]
    #
    # is a list of (unhashable) lists, and other times it returns a list of tupples
    # This needs to be checked in aiida-core before changing the following to a simpler syntax
    obtained_setlist = set()
    for element in builder.instructions.get_dict()[expected_listname]:
        obtained_setlist.add(tuple(element))
    assert obtained_setlist == expected_setlist


@pytest.mark.parametrize('track', [True, False])
@pytest.mark.parametrize('populate_level', [1, 2])
def test_transfer_builder_raise(make_data_source, populate_level, track):
    """Test all raises for the get_transfer_builder utility function."""

    data_source = make_data_source(computer=None, populate_level=populate_level)

    with pytest.raises(ValueError) as execinfo:
        _ = get_transfer_builder(data_source, computer=None, track=track)

    assert 'computer' in str(execinfo.value)
