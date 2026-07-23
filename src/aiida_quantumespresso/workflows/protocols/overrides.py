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

``PwBandsWorkChain`` exposes ``PwBaseWorkChain`` under its ``scf`` and ``bands`` namespaces but
``exclude``s some ports (``clean_workdir`` from both, and ``pw.parent_folder`` from ``bands`` since
the builder wires it from the SCF remote folder). The ``scf``/``bands`` overrides therefore use
*narrowed* ``TypedDict``\\ s (``PwBandsScfOverrides``/``PwBandsBandsOverrides`` and the narrowed
``PwBandsCalculationOverrides``) rather than the full ``PwBaseOverrides``, so autocomplete never
suggests a key the namespace does not expose. The drift guards check each nested namespace against
its own sub-spec so this narrowing cannot silently rot either.
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


class _PwCalculationCommonOverrides(TypedDict, total=False):
    """The ``pw`` (``PwCalculation``) override keys shared by every ``PwBaseWorkChain`` namespace.

    ``code`` and ``structure`` are absent throughout: the builder always sets them from its own
    ``code``/``structure`` arguments, so an override there would be ignored (they are listed in the
    drift guard's ``intentionally_untyped`` set). ``parent_folder`` is *not* here because the
    ``bands`` namespace of ``PwBandsWorkChain`` excludes it; it is added by ``PwCalculationOverrides``.
    """

    parameters: PwParametersOverrides
    pseudos: dict[str, Any]
    metadata: dict[str, Any]
    settings: dict[str, Any]
    parallelization: dict[str, Any]
    monitors: dict[str, Any]
    hubbard_file: Any
    remote_folder: Any
    vdw_table: Any


class PwCalculationOverrides(_PwCalculationCommonOverrides, total=False):
    """Overrides for the ``pw`` namespace of ``PwBaseWorkChain`` (the ``PwCalculation`` inputs)."""

    parent_folder: Any


class PwBandsCalculationOverrides(_PwCalculationCommonOverrides, total=False):
    """Overrides for the ``bands.pw`` namespace of ``PwBandsWorkChain``.

    Identical to ``PwCalculationOverrides`` but without ``parent_folder``: ``PwBandsWorkChain``
    excludes ``pw.parent_folder`` from its ``bands`` namespace (the builder sets it from the SCF
    remote folder), so exposing it here would suggest a key the namespace does not accept.
    """


class PwMetaParameters(TypedDict, total=False):
    """Overrides for the protocol ``meta_parameters`` (consumed by the builder, not an input port)."""

    conv_thr_per_atom: float
    etot_conv_thr_per_atom: float


class _PwBaseCommonOverrides(TypedDict, total=False):
    """The ``PwBaseWorkChain`` override keys shared by the full work chain and its exposed namespaces.

    ``clean_workdir`` and ``pw`` are declared by each concrete subclass: ``PwBandsWorkChain`` excludes
    ``clean_workdir`` from its ``scf``/``bands`` namespaces, and narrows ``pw`` for ``bands``.
    ``meta_parameters`` and ``pseudo_family`` are protocol keys the builder consumes directly (they
    are not input ports); every other key here is a ``PwBaseWorkChain`` input port.
    """

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


class PwBaseOverrides(_PwBaseCommonOverrides, total=False):
    """Overrides accepted by ``PwBaseWorkChain.get_builder_from_protocol``."""

    pw: PwCalculationOverrides
    clean_workdir: bool


class PwBandsScfOverrides(_PwBaseCommonOverrides, total=False):
    """Overrides for the ``scf`` namespace of ``PwBandsWorkChain``.

    Same as ``PwBaseOverrides`` but without ``clean_workdir``, which ``PwBandsWorkChain`` excludes
    from its ``scf`` namespace (cleanup is driven by the parent work chain's own ``clean_workdir``).
    """

    pw: PwCalculationOverrides


class PwBandsBandsOverrides(_PwBaseCommonOverrides, total=False):
    """Overrides for the ``bands`` namespace of ``PwBandsWorkChain``.

    Same as ``PwBandsScfOverrides`` (no ``clean_workdir``) but with ``pw`` further narrowed to
    ``PwBandsCalculationOverrides`` (no ``parent_folder``), matching the ports the ``bands`` namespace
    exposes.
    """

    pw: PwBandsCalculationOverrides


class PwBandsOverrides(TypedDict, total=False):
    """Overrides accepted by ``PwBandsWorkChain.get_builder_from_protocol``.

    ``structure`` is deliberately omitted: the builder always sets it from its own ``structure``
    argument, so an override there would be ignored (see the drift guard's ``intentionally_untyped``).
    """

    scf: PwBandsScfOverrides
    bands: PwBandsBandsOverrides
    clean_workdir: bool
    nbands_factor: float
    bands_kpoints: Any
    bands_kpoints_distance: float
    metadata: dict[str, Any]
