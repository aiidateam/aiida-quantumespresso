# -*- coding: utf-8 -*-
"""A collection of useful pytest fixtures.

* make_tempdir: returns a function to create temporary directories

* makeget_computer: returns a function to create computers (note: these need to be destroyed after the test
  by using 'clear_database' because as the python variables go out of scope, the tempdir used to create the
  computers is deleted).

"""
import pytest

# pylint: disable=redefined-outer-name,bad-option-value,raise-missing-from


# Maybe this could go on aiida-core
@pytest.fixture(scope='function')
def make_tempdir(request):
    """Provide a function to generate new temporary directories.

    :return: The path to the directory
    :rtype: str
    """

    def _make_tempdir(parent_dir=None):
        import tempfile
        import shutil

        try:
            dirpath = tempfile.mkdtemp(dir=parent_dir)
        except FileNotFoundError as exc:
            raise ValueError('The parent_dir provided was never created') from exc

        # After the test function has completed, remove the directory again
        # Design Note: yields would be a simpler way to do this, but can't use it
        # inside a fixture factory without having problems with 'generator' object
        def cleanup():
            shutil.rmtree(dirpath)

        request.addfinalizer(cleanup)

        return dirpath

    return _make_tempdir


# Maybe this could go on aiida-core
# It would be preferrable to do usefixtures but this just works for test, not other
# fixtures (see https://github.com/pytest-dev/pytest/issues/3664).
# We'll have to leave the pylint ignore until then
#@pytest.mark.usefixtures('clear_database')
# pylint: disable=unused-argument
@pytest.fixture(scope='function')
def makeget_computer(clear_database, make_tempdir):
    """Provide a function to generate new computers.

    :return: The computer node
    :rtype: :py:class:`aiida.orm.Computer`
    """

    def _makeget_computer(label='localhost-test', transport_type='local'):
        from aiida.orm import Computer
        from aiida.common.exceptions import NotExistent

        try:
            computer = Computer.objects.get(label=label)

        except NotExistent:
            computer = Computer(
                label=label,
                description=f'{label} computer set up by test manager',
                hostname=label,
                workdir=make_tempdir(),
                transport_type=transport_type,
                scheduler_type='direct'
            )
            computer.store()
            computer.set_minimum_job_poll_interval(0.)

            if transport_type == 'local':
                computer.configure()
            else:
                raise NotImplementedError('Transport `{transport_type}` not implemented')

        return computer

    return _makeget_computer
