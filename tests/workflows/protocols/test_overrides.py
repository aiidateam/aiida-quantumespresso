"""Drift guards for the ``overrides`` ``TypedDict`` definitions.

Each ``TypedDict`` in :mod:`aiida_quantumespresso.workflows.protocols.overrides` mirrors the
``overrides`` mapping consumed by a ``get_builder_from_protocol`` method. These tests assert that the
keys never drift away from the actual work chain input spec, and they do so with *equality* rather
than a one-sided subset:

    set(TheOverrides.__annotations__) == (port_names | builder_consumed) - intentionally_untyped

A subset check only catches the *orphaned* direction — a typed key that no longer maps to any port.
Equality also catches the *incomplete* direction: upstream adds an input port and the ``TypedDict``
silently lacks it, which is exactly the autocomplete/type-check gap this module closes. Because the
overrides recursively merge onto the builder inputs, every input port is legitimate override surface,
so the target set is the full port namespace (plus builder-consumed protocol keys the logic reads
directly without exposing as ports) minus an explicit, justified set of ports the builder always
overwrites from its own arguments. The shared :func:`assert_overrides_match_spec` helper performs the
comparison and reports which direction drifted.
"""

from typing import Callable, Iterable

import pytest
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.workflows.protocols.overrides import (
    PwBandsBandsOverrides,
    PwBandsCalculationOverrides,
    PwBandsProtocolOverrides,
    PwBandsScfOverrides,
    PwBaseProtocolOverrides,
    PwCalculationOverrides,
    PwMetaParameters,
    PwParametersOverrides,
)
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

PwCalculation = CalculationFactory('quantumespresso.pw')

# The complete set of Quantum ESPRESSO ``pw.x`` input namelists that ``PwParametersOverrides``
# enumerates. There is no authoritative in-repo list of *all* ``pw.x`` namelists to tie this to
# (``PwCalculation._automatic_namelists`` only covers the geometry/electronic ones dispatched by
# calculation type; ``FCP``/``RISM`` are QE features aiida-quantumespresso does not reference
# elsewhere). This constant is therefore a second hand-written copy: the equality guard below only
# protects the two declarations against drifting apart, not against QE's real namelist set.
PW_NAMELISTS = {'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'FCP', 'RISM'}

# Protocol keys the ``PwBaseWorkChain`` builder pops directly (``inputs.pop(...)``) without exposing
# them as input ports. The ``scf``/``bands`` overrides of ``PwBandsWorkChain`` are forwarded verbatim
# to that builder, so they carry the same keys.
_BASE_BUILDER_CONSUMED = {'meta_parameters', 'pseudo_family'}


def assert_overrides_match_spec(
    overrides_cls: type,
    port_names: Iterable[str],
    *,
    builder_consumed: Iterable[str] = frozenset(),
    intentionally_untyped: Iterable[str] = frozenset(),
) -> None:
    """Assert the ``TypedDict`` covers exactly the override surface (two-way drift guard).

    :param overrides_cls: the ``TypedDict`` whose ``__annotations__`` keys are the declared overrides.
    :param port_names: the ``spec().inputs`` (or namelist/protocol) key names the type must mirror.
    :param builder_consumed: protocol keys the builder reads directly without exposing as input ports.
    :param intentionally_untyped: ports the builder always overwrites from its own arguments, so an
        override there would be ignored and is deliberately left out of the ``TypedDict``.
    """
    expected = (set(port_names) | set(builder_consumed)) - set(intentionally_untyped)
    actual = set(overrides_cls.__annotations__)
    assert actual == expected, (
        f'missing from {overrides_cls.__name__}: {sorted(expected - actual)}; '
        f'orphaned in {overrides_cls.__name__}: {sorted(actual - expected)}'
    )


# One case per typed namespace, in the argument order of ``assert_overrides_match_spec``: the
# ``TypedDict``, a callable returning the key names it must mirror (called inside the test to keep
# ``spec()`` off the collection path), the builder-consumed keys, and the intentionally untyped ones.
OVERRIDES_CASES = [
    pytest.param(
        PwBaseProtocolOverrides,
        lambda: PwBaseWorkChain.spec().inputs.keys(),
        _BASE_BUILDER_CONSUMED,
        # ``code``/``structure`` land in the nested ``pw`` namespace, so no top-level port is builder-set.
        frozenset(),
        id='base',
    ),
    pytest.param(
        PwCalculationOverrides,
        lambda: PwBaseWorkChain.spec().inputs['pw'].keys(),
        frozenset(),
        {
            'code',  # set unconditionally from the ``code`` argument
            'structure',  # set unconditionally from the ``structure`` argument
        },
        id='base.pw',
    ),
    pytest.param(
        PwBandsProtocolOverrides,
        lambda: PwBandsWorkChain.spec().inputs.keys(),
        frozenset(),
        {
            'structure',  # set unconditionally from the ``structure`` argument
        },
        id='bands',
    ),
    pytest.param(
        PwBandsScfOverrides,
        # The ``scf`` namespace excludes ``clean_workdir``, so the type must not declare it.
        lambda: PwBandsWorkChain.spec().inputs['scf'].keys(),
        _BASE_BUILDER_CONSUMED,
        frozenset(),
        id='bands.scf',
    ),
    pytest.param(
        PwCalculationOverrides,
        lambda: PwBandsWorkChain.spec().inputs['scf']['pw'].keys(),
        frozenset(),
        {
            'code',  # set by the builder (``structure`` is already excluded from the ``scf`` namespace)
        },
        id='bands.scf.pw',
    ),
    pytest.param(
        PwBandsBandsOverrides,
        # The ``bands`` namespace excludes ``clean_workdir``, so the type must not declare it.
        lambda: PwBandsWorkChain.spec().inputs['bands'].keys(),
        _BASE_BUILDER_CONSUMED,
        frozenset(),
        id='bands.bands',
    ),
    pytest.param(
        PwBandsCalculationOverrides,
        # ``PwBandsWorkChain`` excludes ``pw.parent_folder`` (and ``pw.structure``) from its ``bands``
        # namespace, so this type must omit ``parent_folder`` where ``PwCalculationOverrides`` keeps it.
        lambda: PwBandsWorkChain.spec().inputs['bands']['pw'].keys(),
        frozenset(),
        {
            'code',  # set by the builder (``structure``/``parent_folder`` are excluded from the namespace)
        },
        id='bands.bands.pw',
    ),
    pytest.param(
        PwMetaParameters,
        lambda: PwBaseWorkChain.get_protocol_inputs()['meta_parameters'],
        frozenset(),
        # The builder reads every ``meta_parameters`` key, so all of them are override surface.
        frozenset(),
        id='meta_parameters',
    ),
]


@pytest.mark.parametrize('overrides_cls, port_names, builder_consumed, intentionally_untyped', OVERRIDES_CASES)
def test_overrides_match_spec(
    overrides_cls: type,
    port_names: Callable[[], Iterable[str]],
    builder_consumed: Iterable[str],
    intentionally_untyped: Iterable[str],
):
    """The typed overrides keys must equal the override surface of the namespace they mirror."""
    assert_overrides_match_spec(
        overrides_cls,
        port_names(),
        builder_consumed=builder_consumed,
        intentionally_untyped=intentionally_untyped,
    )


def test_pw_parameters_overrides_match_namelists():
    """The ``pw.parameters`` overrides keys must equal the hand-maintained ``PW_NAMELISTS`` copy.

    Both sides are hand-written, so this only guards the two declarations against drifting apart, not
    against Quantum ESPRESSO's real ``pw.x`` namelist set (there is no authoritative in-repo source
    for the latter; see the ``PW_NAMELISTS`` comment). The ``_automatic_namelists`` check below adds a
    partial guard tied to actual plugin code.
    """
    assert_overrides_match_spec(PwParametersOverrides, PW_NAMELISTS)


def test_pw_parameters_overrides_cover_automatic_namelists():
    """Every namelist ``PwCalculation`` dispatches by calculation type must be a typed override key.

    Unlike ``PW_NAMELISTS`` (a second hand-written copy), ``PwCalculation._automatic_namelists`` is
    real plugin code, so this is a genuine one-directional guard: if a calculation type gains a
    namelist that ``PwParametersOverrides`` lacks, this fails. It cannot cover ``FCP``/``RISM``, which
    are not referenced there.
    """
    automatic = set().union(*PwCalculation._automatic_namelists.values())  # noqa: SLF001
    typed = set(PwParametersOverrides.__annotations__)
    assert automatic <= typed, f'namelists dispatched by PwCalculation but not typed: {sorted(automatic - typed)}'
