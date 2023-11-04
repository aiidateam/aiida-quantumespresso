# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
import time

def test_verdi_status(aiida_exec, container_user, timeout):
    time.sleep(timeout)
    output = aiida_exec('verdi status', user=container_user).decode().strip()
    assert 'Connected to RabbitMQ' in output
    assert 'Daemon is running' in output

    # check that we have suppressed the warnings
    assert 'Warning' not in output


def test_computer_setup_success(aiida_exec, container_user, timeout):
    output = aiida_exec('verdi computer test localhost', user=container_user).decode().strip()

    assert "Success" in output
    assert "Failed" not in output

def test_run_real_pw_computation(aiida_exec, container_user, qe_version, sssp_version, timeout):
    import re

    output = aiida_exec("verdi data core.structure import ase /opt/examples/Si.cif", user=container_user).decode().strip()
    
    # Find pk 
    pk = re.search(r"PK = (\d+)", output).group(1)
    
    cmd = f"aiida-quantumespresso calculation launch pw -X pw-{qe_version}@localhost -F SSSP/{sssp_version}/PBE/efficiency -S {pk} -k 1 1 1"
    output = aiida_exec(cmd, user=container_user).decode().strip()

    assert "terminated with state: finished [0]" in output