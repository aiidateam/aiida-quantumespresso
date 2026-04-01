# Migration

## `v5.0.0`

### `PwRelaxWorkChain`

The `PwRelaxWorkChain` has been significantly reworked.
The key changes are:

1. **The `base` namespace is renamed to `base_relax`.**
   Any inputs for the main relaxation loop - previously passed under `base` - must now be passed under `base_relax`.

2. **The final SCF step has been removed.**
   The `base_final_scf` namespace has been removed, and the outputs are taken from the main relax loop instead.

3. **The meta-convergence criterion has changed.**
   The `volume_convergence` input is removed.
   Meta-convergence now checks for Pulay stresses and whether the k-point mesh remains sufficiently dense, rather than comparing cell volumes between iterations.

4. **An optional initial relaxation step has been added.**
   A new `base_init_relax` namespace allows running a first (optional) geometry optimization with looser precision settings, improving robustness for structures far from equilibrium.

```{dropdown} Motivation
:color: info

The final SCF step was redundant: Quantum ESPRESSO already performs one at the end of every `vc-relax`.
The old volume-based meta-convergence criterion always required at least two relaxation runs, even for well-converged structures.
The new criterion — checking for Pulay stresses and k-point mesh density — is both cheaper and more physically meaningful.
The optional initial relaxation step improves robustness for structures far from equilibrium at little extra cost.
```

**Commits**
- ‼️ `PwRelaxWorkChain`: Revisit work chain logic [[`830aef3`](https://github.com/aiidateam/aiida-quantumespresso/commit/830aef35445e8c6d596e54651c9a891ce6b4b77e)]
- 💥 `PwRelaxWorkChain`: remove `volume_convergence` input [[`a31934e`](https://github.com/aiidateam/aiida-quantumespresso/commit/a31934eda0d0fa2867f4fa6c5fe1c10ff97821b6)]

### `PwBandsWorkChain`

The `relax` namespace has been removed from `PwBandsWorkChain`.
The work chain now focuses exclusively on band structure calculations.
If you need to relax the structure first, run `PwRelaxWorkChain` separately and pass the relaxed structure:

```python
relax_builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure)
relax_workchain = engine.submit(relax_builder, wait=True)

# ... wait for workflow to finish

relaxed_structure = relax_workchain.outputs.output_structure
bands_builder = PwBandsWorkChain.get_builder_from_protocol(code, relaxed_structure)
```

```{dropdown} Motivation
:color: info

Combining structure relaxation and band structure calculations in a single workflow conflates two conceptually distinct tasks and creates unnecessary complexity:

- Users who only wanted band structures were forced to understand relaxation parameters, or remove the relaxation namespace when using the protocols.
- Users that do want to relax the structure had to deal with very nested inputs.

Separating the two makes each workflow simpler, more composable, and easier to reason about.

Moreover, if we add the relaxation to `PwBandsWorkChain`, why not to the `PdosWorkChain` as well?
Better to isolate workflows with a certain scope and let the user connect them.
```

**Commits**
- 💥 `PwBandsWorkChain`: remove relaxation [[`6450d7c`](https://github.com/aiidateam/aiida-quantumespresso/commit/6450d7c52c84521ff8d44ad219fff40501e2035b)]

### `parameters`: enforce consistent casing

The `parameters` input of all calculation classes now enforces consistent key casing: namelist keys must be **UPPERCASE** (`CONTROL`, `SYSTEM`, etc.) and parameter keys must be **lowercase**.
Previously, mixed casing was accepted, but logic throughout the codebase would often expect a certain case and break (silently).

```python
# Before — mixed casing was accepted
parameters = Dict({'control': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 60}})

# After — must be consistent
parameters = Dict({'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 60}})
```

If the `Dict` node is not yet stored, a warning is raised and parameters are updated on the fly:

```python
UserWarning: Namelist 'control' should be UPPERCASE. Changed to 'CONTROL'.
```

If the node is stored, a `ValueError` is raised:

```python
ValueError: invalid attribute value Namelist 'control' should be UPPERCASE. Since the Dict node is stored, it cannot be auto-corrected. Please use 'CONTROL' instead.
```

For stored parameter nodes, you will need to create new nodes with the correct casing before passing them to the calculation.

```{dropdown} Motivation
:color: info

FORTRAN namelists are case-insensitive, but Python dictionaries are not.
Allowing mixed casing meant that all internal logic — input file generation, validation, and parameter-based branching — had to defensively handle every possible casing combination.
This made the code verbose, error-prone, and meant that bugs from inconsistent casing could fail silently.
Enforcing a single convention removes this ambiguity and makes queries on stored parameter nodes reliable.
```

**Commits**
- 💥 `parameters`: enforce consistent casing [[`10e10f6`](https://github.com/aiidateam/aiida-quantumespresso/commit/10e10f6c17a053da2f9e8b6629fe611fdb69ebab)]

### Moved packages: XSpectra/XPS and EPW

Support for **XSpectra and XPS** calculations has been removed from `aiida-quantumespresso` and moved to the dedicated [`aiida-qe-xspec`](https://github.com/aiidaplugins/aiida-qe-xspec) package.

Support for **EPW** has also been removed.
The `EpwCalculation` was never fully supported (no parser was implemented), and development has moved to the [`aiida-epw`](https://github.com/aiidaplugins/aiida-epw) package.

```{dropdown} Motivation
:color: info

Hosting XSpectra/XPS and EPW in `aiida-quantumespresso` created maintenance burden for functionality that has a narrow user base and distinct development needs.
For example, having this code in `aiida-quantumespresso` means any change must be backwards-compatible.
Moving them to dedicated packages allows each to evolve independently and keeps `aiida-quantumespresso` focused on the core QE codes.
```

**Commits**
- XSpectra & XPS: Issue Deprecation Warnings [[`c4d7c19`](https://github.com/aiidateam/aiida-quantumespresso/commit/c4d7c190456deca0c26e4fc36304e973fe14a2c2)]
- ‼️ Remove all code related to XSpectra & XPS [[`1e3ea1e`](https://github.com/aiidateam/aiida-quantumespresso/commit/1e3ea1e1ea67dff5d781e86a1afea41e4ec7ac3d)]
- ‼️ Remove EPW-related content [[`63fee4d`](https://github.com/aiidateam/aiida-quantumespresso/commit/63fee4db01cc25ad7ecaec0b958ee661d2ec3de6)]

### Removed legacy parsers

The legacy XML parser modules for Quantum ESPRESSO < v6.2 have been removed.
If you are running calculations with QE < v6.2, you will need to use an older version of `aiida-quantumespresso`.

The `dir_with_bands` parser option for `PwParser` has also been removed, as it is no longer needed.

```{dropdown} Motivation
:color: info

Quantum ESPRESSO v6.2 was released in 2017.
Maintaining parsers for these old XML schemas added complexity and a maintenance burden with little real-world benefit.
```

**Commits**
- 💥 Parsers: Remove legacy code [[`661cda4`](https://github.com/aiidateam/aiida-quantumespresso/commit/661cda4cd911019d13d8a16450c89d0e895bdd0f)]

### Deprecated API

#### `ADDITIONAL_RETRIEVE_LIST`

The `ADDITIONAL_RETRIEVE_LIST` key in the `settings` has been removed, in favor of the `additional_retrieve_list` option for all `CalcJob` classes in `aiida-core`:

```python
inputs['metadata']['options']['additional_retrieve_list'] = ['file1', 'file2']
```

```{dropdown} Motivation
:color: info

`ADDITIONAL_RETRIEVE_LIST` was a plugin-specific workaround predating the `additional_retrieve_list` option that `aiida-core` now provides for any `CalcJob`.
The core option is the proper, standardised way to achieve the same result.
```

#### `PpCalculation` - `keep_plot_file`

The deprecated `keep_plot_file` metadata option of `PpCalculation` has been removed.
Use `keep_data_files` instead:

```python
inputs['metadata']['options']['keep_data_files'] = True
```

```{dropdown} Motivation
:color: info

The `keep_plot_file` name was misleading: the option controls whether intermediate data files are kept, not just plot files.
`keep_data_files` is a more accurate name.
```

#### `NebCalculation` - `first_structure` and `last_structure`

The deprecated `first_structure` and `last_structure` inputs of the `NebCalculation` have been removed.
Pass the initial and final structures via the `images` input as a `TrajectoryData` node instead:

```python
from aiida.orm import TrajectoryData

inputs['images'] = TrajectoryData([initial_structure, final_structure])
```

```{dropdown} Motivation
:color: info

The `images` input as a `TrajectoryData` node is more general — it naturally supports NEB paths with more than two images — and is consistent with how trajectory-like data is handled elsewhere in AiiDA.
```

#### `get_starting_magnetization`

The deprecated `get_starting_magnetization` function has been removed, use `get_magnetization` instead.

```{note}
Both `get_starting_magnetization` and `get_magnetization` were never intended to be public API, but used as utility methods for the (public) `.get_builder_from_protocol()` methods.
However, since we have not strictly defined our public API yet, a deprecation warning was added for `get_starting_magnetization`, and we add a migration guide here.
```

Note that the interface has changed: rather than a `PseudoPotentialFamily` object, `get_magnetization` takes a `z_valences` dictionary mapping element symbols to their valence, for example:

```python
from aiida_quantumespresso.workflows.protocols.utils import get_magnetization

z_valences = {"Si": 4, "O": 6}

magnetization = get_magnetization(structure, z_valences, initial_magnetic_moments)
```

You can still construct the `z_valences` from a `PseudoPotentialFamily` group, for example:

```python
pseudo_family = orm.load_group('SSSP/1.3/PBEsol/efficiency')

z_valences = {
    kind.symbol: pseudo_family.get_pseudo(element=kind.symbol).z_valence
    for kind in structure.kinds
}
```

```{dropdown} Motivation
:color: info

`get_starting_magnetization` only handled collinear calculations.
`get_magnetization` was introduced to also support non-collinear and spin-orbit calculations via the `spin_type` parameter, returning `angle1` and `angle2` in addition to `starting_magnetization`.
A subsequent refactor also replaced the `PseudoPotentialFamily` argument with `z_valences`, a plain dictionary, making the function usable with custom pseudopotentials.
```

**Commits**
- 🗑️ Deprecate `ADDITIONAL_RETRIEVE_LIST` in `settings` [[`5d51932`](https://github.com/aiidateam/aiida-quantumespresso/commit/5d51932b6fb1e964e0753ccce32d40a35fb8ba89)]
- 💥 Calculations: Remove `ADDITIONAL_RETRIEVE_LIST` setting [[`23072e1`](https://github.com/aiidateam/aiida-quantumespresso/commit/23072e1)]
- 💥 `PpCalculation`: Remove `keep_plot_file` option [[`25e7dcc`](https://github.com/aiidateam/aiida-quantumespresso/commit/25e7dcc)]
- 👌 `NebCalculation`: Use `images` instead of `first/last_structure` [[`b0b63c1`](https://github.com/aiidateam/aiida-quantumespresso/commit/b0b63c1c2493949a8ffe95fe59ae2da98ffe6b62)]
- 💥 `NebCalculation`: Remove `first_structure`/`last_structure` inputs [[`c4da42a`](https://github.com/aiidateam/aiida-quantumespresso/commit/c4da42a)]
- ✨ Protocols: Add support for non-collinear calculations [[`5586c3f`](https://github.com/aiidateam/aiida-quantumespresso/commit/5586c3fee9bd4118c7a61b52b913540db356832e)]
- 👌 `get_magnetization`: Replace `pseudo_family` input by `z_valences` [[`16fbada`](https://github.com/aiidateam/aiida-quantumespresso/commit/16fbadafe5bd240004f173c506941e4866d08320)]
- 💥 Protocols: Remove `get_starting_magnetization` function [[`673a970`](https://github.com/aiidateam/aiida-quantumespresso/commit/673a970)]
