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

from typing import Iterable

from aiida_quantumespresso.workflows.protocols.overrides import (
    PwBandsOverrides,
    PwBaseOverrides,
    PwCalculationOverrides,
    PwMetaParameters,
    PwParametersOverrides,
)
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

# The complete set of Quantum ESPRESSO ``pw.x`` input namelists. Any of these may be merged into the
# ``pw.parameters`` ``Dict`` through the overrides, so ``PwParametersOverrides`` enumerates them all.
PW_NAMELISTS = {'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'FCP', 'RISM'}


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


def test_pw_parameters_overrides_match_namelists():
    """The ``pw.parameters`` overrides keys must be exactly the Quantum ESPRESSO ``pw.x`` namelists."""
    assert_overrides_match_spec(PwParametersOverrides, PW_NAMELISTS)


def test_pw_meta_parameters_match_protocol():
    """The ``meta_parameters`` overrides keys must match the protocol ``meta_parameters`` schema."""
    protocol_meta = PwBaseWorkChain.get_protocol_inputs()['meta_parameters']
    # No exclusions: the builder reads every ``meta_parameters`` key, so all are override surface.
    assert_overrides_match_spec(PwMetaParameters, protocol_meta)


def test_pw_calculation_overrides_match_spec():
    """The ``pw`` overrides keys must equal the ``pw`` input namespace ports (bar builder-set ones)."""
    ports = PwBaseWorkChain.spec().inputs['pw'].keys()
    # Ports the builder always sets from its own arguments, so an override would be ignored.
    intentionally_untyped = {
        'code',  # set unconditionally from the ``code`` argument
        'structure',  # set unconditionally from the ``structure`` argument
    }
    assert_overrides_match_spec(PwCalculationOverrides, ports, intentionally_untyped=intentionally_untyped)


def test_pw_base_overrides_match_spec():
    """The ``PwBaseWorkChain`` overrides keys must equal spec ports plus builder-consumed keys."""
    ports = PwBaseWorkChain.spec().inputs.keys()
    # Protocol keys the builder pops directly (``inputs.pop(...)``) without exposing as input ports.
    builder_consumed = {'meta_parameters', 'pseudo_family'}
    # Nothing at this level is set from the builder's own arguments (``code``/``structure`` land in
    # the nested ``pw`` namespace), so no port is excluded.
    assert_overrides_match_spec(PwBaseOverrides, ports, builder_consumed=builder_consumed)


def test_pw_bands_overrides_match_spec():
    """The ``PwBandsWorkChain`` overrides keys must equal spec ports (nested ``scf``/``bands`` recurse)."""
    ports = PwBandsWorkChain.spec().inputs.keys()
    # Ports the builder always sets from its own arguments, so an override would be ignored.
    intentionally_untyped = {
        'structure',  # set unconditionally from the ``structure`` argument
    }
    assert_overrides_match_spec(PwBandsOverrides, ports, intentionally_untyped=intentionally_untyped)
