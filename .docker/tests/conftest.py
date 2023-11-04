# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring, redefined-outer-name
import json
import re
from pathlib import Path

import pytest


@pytest.fixture(scope='session', params=['aiida-quantumespresso'])
def variant(request):
    return request.param


@pytest.fixture(scope='session')
def docker_compose_file(pytestconfig, variant):  # pylint: disable=unused-argument
    return f'docker-compose.{variant}.yml'


@pytest.fixture(scope='session')
def docker_compose(docker_services):
    # pylint: disable=protected-access
    return docker_services._docker_compose


@pytest.fixture
def timeout():
    """Container and service startup timeout"""
    return 60


@pytest.fixture
def container_user():
    return 'aiida'


@pytest.fixture
def aiida_exec(docker_compose):

    def execute(command, user=None, **kwargs):
        if user:
            command = f'exec -T --user={user} aiida {command}'
        else:
            command = f'exec -T aiida {command}'
        return docker_compose.execute(command, **kwargs)

    return execute

@pytest.fixture
def qe_version(aiida_exec):
    info = json.loads(
        aiida_exec(
            "mamba list -n base --json --full-name qe"
        ).decode()
    )[0]
    return info['version']

@pytest.fixture
def sssp_version(aiida_exec):
    output = aiida_exec("aiida-pseudo list").decode().strip()
    return re.search(r"SSSP/(\d+\.\d+)/PBE/efficiency", output).group(1)
