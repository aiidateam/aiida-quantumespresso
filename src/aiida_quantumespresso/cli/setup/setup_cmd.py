# -*- coding: utf-8 -*-
"""Command line interface commands for setting up `code`s for Quantum ESPRESSO executables."""
# -*- coding: utf-8 -*-
import pathlib

from aiida import orm
from aiida.cmdline.params import arguments
from aiida.cmdline.params.options.commands import code as options_code
from aiida.cmdline.utils import echo
from aiida.common.exceptions import NotExistent
import click


# add option to specify executables, otherwise try to set all that can be found in the directory
# documentation: this is equivalen to: list of executables e.g. pw.x dos. projwfc.x ...

# suffix: try to use curly brackets and f-strings {}_, show example in docstrings

# version in the name
# try to run verdi code test and parse version number
# --test option to run the test of code
# maybe, automatic test does only work if we don't submit
# just use transport to run code

@arguments.COMPUTER()
@click.argument('executables', nargs=-1, required=True)
@click.option('--directory', '-d', help='Absolute path to directory where the executbale(s) is (are) located.')
@click.option('--label_suffix', '-l', help='Suffix to append to the label of the code instance.', default='')
@options_code.PREPEND_TEXT()
@options_code.APPEND_TEXT()
def create_code(executables, computer, directory, label_suffix, **kwargs):
    """Create a `code` instance for a Quantum ESPRESSO executable."""
    user = orm.User.collection.get_default()

    try:
        authinfo = computer.get_authinfo(user)
    except NotExistent:
        echo.echo_critical(f'Computer<{computer.label}> is not yet configured for user<{user.email}>')

    if not authinfo.enabled:
        echo.echo_warning(f'Computer<{computer.label}> is disabled for user<{user.email}>')
        click.confirm('Do you really want to test it?', abort=True)

    transport = authinfo.get_transport()

    prepend_text = kwargs.pop('prepend_text', None)
    executable_paths = []
    with transport:
        for executable in executables:
            if not directory:
                if prepend_text:
                    transport.exec_command_wait(prepend_text)

                # Add prepend text to the command? Sometimes one might need to use modules before finding the executable
                which_ret_val, exec_path, which_stderr = transport.exec_command_wait(f'which {executable}')

                if which_ret_val != 0:
                    raise ValueError(
                        f'Failed to determine the path of the executable<{executable}> on '
                        f'computer<{computer.label}>: {which_stderr}'
                    )
            else:
                directory = pathlib.PurePosixPath(directory)
                if not directory.is_absolute():
                    raise ValueError(f'Directory<{directory}> is not an absolute path.')
                if not transport.path_exists(directory):
                    raise ValueError(f'Directory<{directory}> does not exist on computer<{computer.label}>')
                exec_path = directory / executable
            # also check that executable is executable and exists
            executable_paths.append(str(exec_path))

    for executable, exec_path in zip(executables, executable_paths):
        label = f'{executable.replace(".x", "")}_{label_suffix}' if label_suffix else executable.replace(".x", "")
        code = orm.InstalledCode(
            label=label,
            computer=computer,
            filepath_executable=exec_path,
            default_calc_job_plugin=f'quantumespresso.{executable.split(".")[0]}',
            prepend_text=prepend_text
        )
        code.store()
        echo.echo_success(f'Code<{code.label}> for {executable} created with pk<{code.pk}> on computer<{computer.label}>')

# add option --no-connect, makes -d mandatory
# In case you can not connect

# Run test, start code
# parse general information of the code e.g. simply run `pw.x` and append to extras
# extras can also be used to validate if things change

# aiida-quantumespresso code test
    # run on login node
    # --via-submit --> check report when it's done
    # set extras on the code in the end
    
# set extras, option -- don't set extras
    # md5 of the file
    # version of the code
    # compiled features? e.g. scalapack, hdf5, ...
    # If md5 changes, fails --> optional force extras update
    # Fails if others change as well
    
# test multiple codes
# use similar template from setup to run multiple codes that match the same pattern
# unless disabled, will also test code version and add it to the extras
# --provide-relabel-command/relabel-pattern=
# final report will show how to relabel codes based on the found information
# make it easily accessible to copy and paste

# relabel in interactive mode --relabel, -n option
# What to do if code with given label already exists?
# report existing codes and incrementing index, provide a command
# that could be run to add this

# check which spearator to use, dash might be confusing

# @cmd_setup.command('test')
# def test_code():
#     """Test a `code` instance for a Quantum ESPRESSO executable."""
#     from ase.build import bulk

#     # Create a silicon crystal structure
#     structure = bulk('Si')
#     structure = orm.StructureData(ase=structure)

#     # Create a Quantum ESPRESSO calculation
#     calculation = orm.CalcJobNode()
