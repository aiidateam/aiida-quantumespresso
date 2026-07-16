"""Typed dictionaries describing the ``overrides`` accepted by ``get_builder_from_protocol``.

These ``TypedDict`` definitions mirror the nested ``overrides`` mapping that each work chain's
``get_builder_from_protocol`` recursively merges onto the protocol defaults. They exist purely to
provide editor autocompletion and static type checking for callers; they carry no runtime behaviour
and the builder logic remains unchanged whether or not they are used.

To keep them from silently drifting away from the work chain input specs, the keys of each
``TypedDict`` are asserted to be *equal* to the corresponding ``spec().inputs`` port names (plus a
small set of "builder-consumed" keys that the protocol logic reads directly without exposing them as
input ports, minus an explicit set of ports the builder always overwrites from its own arguments) in
``tests/workflows/protocols/test_overrides.py``. Equality (rather than a one-sided subset) means the
guard fires in both directions: an orphaned typed key *and* a newly added upstream port that the
``TypedDict`` is missing both break the test.
"""

from __future__ import annotations

from typing import Any, TypedDict


class PwParametersOverrides(TypedDict, total=False):
    """Overrides for the ``pw.parameters`` Quantum ESPRESSO namelists.

    The keys are the QE ``pw.x`` namelist names; the values are free-form mappings of namelist
    keyword to value. This enumerates the full ``pw.x`` namelist set so that any namelist a caller
    may legitimately merge in through the ``overrides`` is autocompleted.
    """

    CONTROL: dict[str, Any]
    SYSTEM: dict[str, Any]
    ELECTRONS: dict[str, Any]
    IONS: dict[str, Any]
    CELL: dict[str, Any]
    FCP: dict[str, Any]
    RISM: dict[str, Any]


class PwCalculationOverrides(TypedDict, total=False):
    """Overrides for the ``pw`` namespace of ``PwBaseWorkChain`` (the ``PwCalculation`` inputs).

    ``code`` and ``structure`` are deliberately omitted: the builder always sets them from its own
    ``code``/``structure`` arguments, so an override there would be ignored (see the drift-guard's
    ``INTENTIONALLY_UNTYPED`` set).
    """

    parameters: PwParametersOverrides
    pseudos: dict[str, Any]
    metadata: dict[str, Any]
    settings: dict[str, Any]
    parallelization: dict[str, Any]
    monitors: dict[str, Any]
    hubbard_file: Any
    parent_folder: Any
    remote_folder: Any
    vdw_table: Any


class PwMetaParameters(TypedDict, total=False):
    """Overrides for the protocol ``meta_parameters`` (consumed by the builder, not an input port)."""

    conv_thr_per_atom: float
    etot_conv_thr_per_atom: float


class PwBaseOverrides(TypedDict, total=False):
    """Overrides accepted by ``PwBaseWorkChain.get_builder_from_protocol``.

    ``meta_parameters`` and ``pseudo_family`` are protocol keys the builder consumes directly (they
    are not input ports); every other key is a ``PwBaseWorkChain`` input port.
    """

    pw: PwCalculationOverrides
    clean_workdir: bool
    kpoints: Any
    kpoints_distance: float
    kpoints_force_parity: bool
    max_iterations: int
    metadata: dict[str, Any]
    handler_overrides: dict[str, Any]
    on_unhandled_failure: str
    pause_on_max_iterations: bool
    meta_parameters: PwMetaParameters
    pseudo_family: str


class PwBandsOverrides(TypedDict, total=False):
    """Overrides accepted by ``PwBandsWorkChain.get_builder_from_protocol``.

    ``structure`` is deliberately omitted: the builder always sets it from its own ``structure``
    argument, so an override there would be ignored (see the drift-guard's ``INTENTIONALLY_UNTYPED``).
    """

    scf: PwBaseOverrides
    bands: PwBaseOverrides
    clean_workdir: bool
    nbands_factor: float
    bands_kpoints: Any
    bands_kpoints_distance: float
    metadata: dict[str, Any]
