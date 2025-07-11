## v4.11.0

Minor release to add support for Quantum ESPRESSO v7.4.1 and set a lower limit for `aiida-pseudo`, due to a change in the SSSP download URL on the Materials Cloud.
A little bug was also fixed for the `pp.x` parser.
The `ADDITIONAL_RETRIEVE_LIST` setting was also deprecated in favor of using the `additional_retrieve_list` option that is now available for any `CalcJob`, for example:

```python
from aiida_quantumespresso.calculations.pw import PwCalculation

builder = PwCalculation.get_builder()
builder.metadata.options['additional_retrieve_list'] = ['filename']
```

### ✨ New features

* `PwParser`: Add XML schema for Quantum ESPRESSO v7.4.1 [[4a6cdc2](https://github.com/aiidateam/aiida-quantumespresso/commit/4a6cdc22ec7f868a497abfd384bfeed39dfee3dc)]

### 🗑️ Deprecations

* Deprecate `ADDITIONAL_RETRIEVE_LIST` in `settings` [[5d51932](https://github.com/aiidateam/aiida-quantumespresso/commit/5d51932b6fb1e964e0753ccce32d40a35fb8ba89)]

### 👌 Improvements

* Fix error message for `ERROR_CHARGE_IS_WRONG` [[ab3318b](https://github.com/aiidateam/aiida-quantumespresso/commit/ab3318bf96aba84a453283406e32f3dc7ba7ac85)]

### 🐛 Bug fixes

* `PpParser`: Fix logic for LDOS Files [[1385c49](https://github.com/aiidateam/aiida-quantumespresso/commit/1385c49e91d3c3c4f7d1dd48d72ce30c044b25bd)]

### 📚 Documentation

* Add tutorial for Hubbard corrections  [[f18961a](https://github.com/aiidateam/aiida-quantumespresso/commit/f18961a02757ff093e41ebe800d55e766b303680)]
* `PwBandsWorkChain`: Fix typo in docstring [[833631c](https://github.com/aiidateam/aiida-quantumespresso/commit/833631c369544be61afb6f76a33997c48b99cd51)]

### ⬆️ Update dependencies

* Bump aiida-pseudo to v1.7  [[bd98428](https://github.com/aiidateam/aiida-quantumespresso/commit/bd984289297e8d353c498976434b980a9d853ab6)]
* Add upper limit to `aiida-pseudo` dependency [[904cbee](https://github.com/aiidateam/aiida-quantumespresso/commit/904cbee313397812b6e7eb7fdb4d0634e483366e)]

### 🧪 Tests

* Clean up `pp.x` 3D parsing failure tests [[dbf9651](https://github.com/aiidateam/aiida-quantumespresso/commit/dbf9651d969e7b90c6e8041197a443e38f39056a)]

## v4.10.0

### 👌 Improvements

* Protocols: Update k-point sampling protocols [[17f0fd8](https://github.com/aiidateam/aiida-quantumespresso/commit/17f0fd893470abebf5d84015be4586d6345252fd)]
* XSpectra & XPS: Issue Depreciation Warnings [[c4d7c19](https://github.com/aiidateam/aiida-quantumespresso/commit/c4d7c190456deca0c26e4fc36304e973fe14a2c2)]

### 🐛 Bug fixes

* `PdosWorkChain`: Enable energy range around the Fermi level [[435fc90](https://github.com/aiidateam/aiida-quantumespresso/commit/435fc9008d012c9a0d81d8fad95624dde7b670cb)]

### 🔧 Maintenance

* CD & nightly: update Github actions versions [[f356d3a](https://github.com/aiidateam/aiida-quantumespresso/commit/f356d3a0c347ade636e75ce063ac0b9c69e5cda6)]
* Update GitHub actions versions [[91c7c30](https://github.com/aiidateam/aiida-quantumespresso/commit/91c7c30bc5b05a5ee85cf0cea9fbed784305ff59)]

### 🧪 Tests

* Create two pseudo families in `sssp` fixture [[043664e](https://github.com/aiidateam/aiida-quantumespresso/commit/043664e2ab607fe922a1eba2b0ee104f50a9616a)]


## v4.9.0

### 👌 Improvements

* `PdosWorkChain`: Update protocols for improved accuracy [[8a13266](https://github.com/aiidateam/aiida-quantumespresso/commit/8a13266be795e31b6ce6c255612596f1377d4478)]

### 🐛 Bug fixes

* `ProjwfcParser`: Fix DOS parsing for spin-polarised case [[bcd5508](https://github.com/aiidateam/aiida-quantumespresso/commit/bcd550847b0bf8142bb43425a1480edfd4db2965)]

### 🧪 Tests

* Allow glob patterns in `retrieve_temporary` input [[af3739c](https://github.com/aiidateam/aiida-quantumespresso/commit/af3739c908e59ae3f11dd78f232703502ef1f0df)]

### ♻️ Refactor

* Refactor `projwfc.x` parser [[29f9c98](https://github.com/aiidateam/aiida-quantumespresso/commit/29f9c9846b9ee63e25f7e28664974228c45ec97c)]

## v4.8.0

### ✨ New features

* `Hubbard`: add useful utility functions  [[41f196b](https://github.com/aiidateam/aiida-quantumespresso/commit/41f196b2b46022c179bdce0dcd44ab338d956551)]
* `PwParser`: Add XML schema for Quantum ESPRESSO v7.4 [[75808a9](https://github.com/aiidateam/aiida-quantumespresso/commit/75808a918075b7ef0ac3b82557b5325330ba0fdc)]

### 👌 Improvements

* `PhParser`: allow for pattern initialization  [[74a18bc](https://github.com/aiidateam/aiida-quantumespresso/commit/74a18bc59beafd62082b24f809bfba29a402489d)]
* `PwBaseWorkChain`: Always do full restart for `ERROR_OUT_OF_WALLTIME` [[fcb8da9](https://github.com/aiidateam/aiida-quantumespresso/commit/fcb8da9a165fc4f6e94c29a6e4f1d29305da4548)]

### 📚 Documentation

* `README.md`: Correct `aiida-core` compatibility [[1c8a804](https://github.com/aiidateam/aiida-quantumespresso/commit/1c8a8049c005e791300b7d2ca7c4beae32257dc9)]

### 🔧 Maintenance

* Devops: Add explicit sphinx.configuration key to RTD conf  [[733c43f](https://github.com/aiidateam/aiida-quantumespresso/commit/733c43ff6b690c19eacb108154038ca93082d145)]
* Update Python support: drop v3.8 and add v3.12 [[172540b](https://github.com/aiidateam/aiida-quantumespresso/commit/172540be5078e77d58b6a8cc1b28cbae2c6e70b9)]

## v4.7.0

### ✨ New features

* Add `nbands_factor` logic into PdosWorkChain  [[1020b02](https://github.com/aiidateam/aiida-quantumespresso/commit/1020b02c76bd3ae9783087bdf5f796380a7fdf3b)]
* ✨ `PwParser`: Add the XML schema for Quantum ESPRESSO v7.3.1 [[57e7463](https://github.com/aiidateam/aiida-quantumespresso/commit/57e7463c5727775d6a0470a41d1aca0ec4083b9a)]
* `XspectraCrystalWorkChain`: Enable Symmetry Data Inputs  [[b79189d](https://github.com/aiidateam/aiida-quantumespresso/commit/b79189d7ce4756e846ab39c567ba4681474741ed)]
* Add calcjob, parser and base workchain plugin for `bands.x`  [[651fd01](https://github.com/aiidateam/aiida-quantumespresso/commit/651fd0142a965ca1b03cc52f0f2f8d960936a1cd)]

### 👌 Improvements

* `PpCalculation`: Make parsing of output files optional  [[bc0d815](https://github.com/aiidateam/aiida-quantumespresso/commit/bc0d8156f3f206b76e15f0f0c0742d8b579b4722)]

### 🐛 Bug fixes

* CLI: Fix bug in `aiida-quantumespresso workflow launc pw-base`  [[ea76d9b](https://github.com/aiidateam/aiida-quantumespresso/commit/ea76d9b37f78315bbf93f93fa56460c7dfe0652a)]

### 📚 Documentation

* Docs: Fix build by pinning `sphinx-autoapi~=3.0.0`  [[91c3e1d](https://github.com/aiidateam/aiida-quantumespresso/commit/91c3e1d35939491663a697d201dcccdf90c076c6)]

### ♻️ Refactor

* `get_xspectra_structures`: Refactor and Improve Code  [[210c40b](https://github.com/aiidateam/aiida-quantumespresso/commit/210c40bbc3445f55155bbb855d320afa00fa347e)]

## v4.6.0

This minor release provides several improvements and bug fixes, mostly related to the `HubbardStructureData` and XPS/XAS calculations.

Since there were no changes in the schema for Quantum ESPRESSO v7.3, versions 4.3 and above of the plugin package should now also fully support the new Quantum ESPRESSO release.

### 👌 Improvements

* XAS: Enable Correct Parsing of Hubbard and Magnetic Data [[f439504](https://github.com/aiidateam/aiida-quantumespresso/commit/f4395048dfb9c74b97a5d38eff99029449816dc0)]
* `PhCalculation`: add symmetry related exit codes [[5a6529f](https://github.com/aiidateam/aiida-quantumespresso/commit/5a6529f46fa3519f006527c02db3a065f6e5ebaa)]
* `seekpath_structure_analysis`: `HubbardStructureData` compatibility [[9cb1cfa](https://github.com/aiidateam/aiida-quantumespresso/commit/9cb1cfa8a70d19af7aaa1b624cf17c8babe93f41)]
* `HubbardStructureData`: add compatibility for <3D structures [[d645069](https://github.com/aiidateam/aiida-quantumespresso/commit/d645069d1d6ac40d6a23e4bbf49b01b51f5bb33c)]
* `HubbardStructureData`: Add validation on site indices [[960a371](https://github.com/aiidateam/aiida-quantumespresso/commit/960a371d2e7f327a3f5bedb86e21dcb5100c03d1)]

### 🐛 Bug fixes

* `PwBandsWorkChain`: Respect `bands_kpoints` in overrides  [[ae7d248](https://github.com/aiidateam/aiida-quantumespresso/commit/ae7d2484a084ca55dcef4c1ca332b2fb0b478fad)]
* CLI: Fix import of `StructureData` from QE input file  [[3440623](https://github.com/aiidateam/aiida-quantumespresso/commit/34406230023c3c3715264a7de180a576fd7def48)]
* Fix inputs for molecules in the XPS calculation  [[b7f17cf](https://github.com/aiidateam/aiida-quantumespresso/commit/b7f17cfaf71e5b07089426a4fbe7ae2ca5317523)]

### 🧪 Tests

* CLI: Add test for data structure import [[5c3c301](https://github.com/aiidateam/aiida-quantumespresso/commit/5c3c3012c36577cc0ca6a4374fd97f2ce3177ebe)]


## v4.5.1

This patch release fixes some issues with the changes introduced in [b9c7517](https://github.com/aiidateam/aiida-quantumespresso/commit/b9c7517744e645a93d4afc9b1999881fc39a0e46) and released in v4.5.0.
The new approach for setting the q-points was unfortunately broken, which is now fixed in [c353cc2](https://github.com/aiidateam/aiida-quantumespresso/commit/c353cc2e4104352ef9b5490adb53a60da47f293d).
Moreover, the validation that was added to the top-level inputs of the `PhBaseWorkChain` requires the user to specify either the `qpoints` or `qpoints_distance` input.
This means that work chains which wrap the `PhBaseWorkChain` but provide the q-points on the fly will have to disable this validation by setting the corresponding validator to `None` in the input spec.
To avoid this, we have the validator check if one of the `qpoints` or `qpoints_distance` ports are still present in the port namespace.
If not, the validation is skipped.

Higher-level work chains that wrap the `PhBaseWorkChain` can then simply exclude these ports when exposing the inputs:

```python
    class WrapPhBaseWorkChain(WorkChain):
        """Example work chain that wraps a ``PhBaseWorkChain`` excluding q-points inputs."""

        @classmethod
        def define(cls, spec):
            super().define(spec)
            spec.expose_inputs(PhBaseWorkChain, exclude=('qpoints', 'qpoints_distance'))
```

### 👌 Improvements

* `PhBaseWorkChain`: skip q-points validation if ports are excluded [[32536e8](https://github.com/aiidateam/aiida-quantumespresso/commit/32536e85abd6de30cd8f9a07124996ce8cd0760a)]

### 🐛 Bug fixes

* `PhBaseWorkChain`: fix `set_qpoints` step [[c353cc2](https://github.com/aiidateam/aiida-quantumespresso/commit/c353cc2e4104352ef9b5490adb53a60da47f293d)]

### 📚 Documentation

* `CHANGELOG.md`: improve release notes for `v4.5.0` [[b659625](https://github.com/aiidateam/aiida-quantumespresso/commit/b65962565f300fbda643ea43d18d01329c4a85ff)]

## v4.5.0

Besides several bug fixes and documentation improvements, this minor release introduces some changes to the `PhBaseWorkChain` and how it is used.
Most importantly, q-points are no longer defined directly on the `PhCalculation` in the `ph` namespace, but as a top-level input on the `PhBaseWorkChain`, either as a `KpointsData` node (see [b9c7517](https://github.com/aiidateam/aiida-quantumespresso/commit/b9c7517744e645a93d4afc9b1999881fc39a0e46)):

```python
ph_base_builder = PhBaseWorkChain.get_builder()

qpoints = orm.KpointsData()
qpoints.set_kpoints_mesh([2, 2, 2])

ph_base_builder.qpoints = qpoints
```

Or a `qpoints_distance`, which defines a linear density that can be used in combination with the structure to generate a mesh on the fly:

```python
ph_base_builder = PhBaseWorkChain.get_builder()

ph_base_builder.qpoints_distance = orm.Float(0.3)
```

The protocols are also updated to use the new `qpoints_distance` input, with some reasonable default for the various options.

| ⚠️ **Important**: While the current defaults are reasonable, they by no means represent a rigorously tested protocol that guarantees a certain level of convergence. Be sure to run your own tests! |
|---|

The `get_builder_from_protocol()` method of the `PhBaseWorkChain` now also recognises the `electronic_type` input argument.
For example, when running an insulator, one can obtain a fully populated builder via:

```python
from aiida_quantumespresso.common.types import ElectronicType

builder = PhBaseWorkChain.get_builder_from_protocol(
    code=ph_code,
    parent_folder=pw_remote_folder,
    electronic_type=ElectronicType.INSULATOR
)
```

which will set `INPUTPH.epsil` to `True` in the `ph.x` input file.

This release also bumps the default SSSP version in the protocols up to v1.3 (still using the PBEsol functional).
In case you have not installed this newer version, you can do so with `aiida-pseudo` (version `>=1.1.0`):

```console
aiida-pseudo install sssp -v 1.3 -x PBEsol
```

### ‼️ Breaking changes

- `PhBaseWorkChain`: allow generation of q-point mesh via `qpoints_distance` [[b9c7517]](https://github.com/aiidateam/aiida-quantumespresso/commit/b9c7517744e645a93d4afc9b1999881fc39a0e46)

### ✨ New features

- `PwBaseWorkChain`: new handler for BFGS history failure [[0224f8a]](https://github.com/aiidateam/aiida-quantumespresso/commit/0224f8a4cf8122a916e2b2c5be11f3a8a811f740)
- Support calculating the XPS spectra of the atoms specific by indices [[fc1a940]](https://github.com/aiidateam/aiida-quantumespresso/commit/fc1a940d4a60f22b42ec0a069a06436c0c9ae0f5)

### 👌 Improvements

- Improve `PhBaseWorkChain` overrides/protocol [[39287e0]](https://github.com/aiidateam/aiida-quantumespresso/commit/39287e03cb6bbf1915662685a5c441e9c7c36030)
- `PwBaseWorkChain`: Remove disabling of resource validation [[d4e6681]](https://github.com/aiidateam/aiida-quantumespresso/commit/d4e668195d369e360bfb1f06611049a940640843)
- Protocols: Bump default SSSP version to 1.3 [[49d503d]](https://github.com/aiidateam/aiida-quantumespresso/commit/49d503d8b2a0c09dd2b38fecb73b28c82e930822)

### 🐛 Bug fixes

- `PdosWorkChain`: Fix constrained magnetization case [[a68e1e1]](https://github.com/aiidateam/aiida-quantumespresso/commit/a68e1e15c11f6ad4461921145c648d75c49ff26c)
- Fix missing `max_iterations` in overrides [[9061ea5]](https://github.com/aiidateam/aiida-quantumespresso/commit/9061ea5df65e4fe95f77309d9abb5a3c7f64bb9f)
- `OpenGridCalculation`: Add the `output_parameters` output to the spec [[5f0e095]](https://github.com/aiidateam/aiida-quantumespresso/commit/5f0e095647d6529c295002aa15e48ff647111ab7)
- `Q2rCalculation`: Add the `output_parameters` output to the spec [[7a303f9]](https://github.com/aiidateam/aiida-quantumespresso/commit/7a303f9e10ff6dd9eae652b4f2d8ad2c482022d6)
- `PwBaseWorkChain`: Pop `starting_magnetization` if `tot_magnetization` is defined [[2adf033]](https://github.com/aiidateam/aiida-quantumespresso/commit/2adf0335291b452d41ee9282ef13ed21f47cebc8)
- `PwCalculation`: Fix calling input validation of base class [[17e173f]](https://github.com/aiidateam/aiida-quantumespresso/commit/17e173f11c75142755bc9f3c9a71160d5aba778c)

### 📚 Documentation

- Remove `aiida.manage.configuration.load_documentation_profile` [[f1d547c]](https://github.com/aiidateam/aiida-quantumespresso/commit/f1d547c28b35241a53703b0790cf5dac70455060)
- Address warning from `pydata-sphinx-theme` [[74bbaa2]](https://github.com/aiidateam/aiida-quantumespresso/commit/74bbaa22b383b3323fcc3d41ad5b82fa89895c92)
- Docs: Fix broken build by updating `sphinx-autoapi~=3.0` [[80e550e]](https://github.com/aiidateam/aiida-quantumespresso/commit/80e550e9e4d831b620f61fb93c88fcc0778a467d)
- Docs: Update QE compatibility matrix in README.md [[5db3b28]](https://github.com/aiidateam/aiida-quantumespresso/commit/5db3b28ca5e067e63e59fdfdf6be8362efc7d223)

### 🔧 Maintenance

- Address various deprecation warnings from `aiida-core` [[f133b9a]](https://github.com/aiidateam/aiida-quantumespresso/commit/f133b9ab8c87c122e0edff2d27bad54d5d834681)

### ⬆️ Update dependencies
- Dependencies: Update `pydantic~=2.4` [[740e0be]](https://github.com/aiidateam/aiida-quantumespresso/commit/740e0bec0e68b3229367b9f10b181a925616c08b)
- Dependencies: Update `xmlschema~=2.0` [[bec6dd6]](https://github.com/aiidateam/aiida-quantumespresso/commit/bec6dd6b56b4cd3bbed3f3ab8fb97c7f7bdc0214)

### 🧪 Tests

- Revert change to `fixture_code` [[e9ce7a0]](https://github.com/aiidateam/aiida-quantumespresso/commit/e9ce7a069cf1011f829143563730750a9b9fc637)


## v4.4.0

This minor release mainly includes documentation changes, but also a small breaking change in [[a389629](https://github.com/aiidateam/aiida-quantumespresso/commit/a389629387b74805ffe2f4d6515ac05b8f62b4d5)].
Work chains that wrap the `PwBaseWorkChain` for an `nscf` or `bands` calculation and add the `parent_folder` during runtime would have to adapt the `validator` for the `pw` name space as follows to avoid warnings:

```python
    spec.inputs['nscf']['pw'].validator = PwCalculation.validate_inputs_base
```

Now, the validator that checks for the presence of the `parent_folder` in the inputs is adapted so it checks if the  `parent_folder` port is in the name space.
If not, the validation is skipped.
So work chains that wrap the `PwBaseWorkChain` to run an `nscf` or `bands` calculation can simply exclude the `parent_folder` port in the `pw` name space:

```python
    spec.expose_inputs(
        PwBaseWorkChain,
        namespace='nscf',
        exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
        namespace_options={
            'help': 'Inputs for the `PwBaseWorkChain` of the `nscf` calculation.',
            'validator': validate_nscf
        }
    )
```

A new calculation function `create_magnetic_configuration` is added, which can be used to create a new `StructureData` with the required kinds for a specific magnetic configuration.
For example, for an HCP structure with two Co sites:

```python
In [1]: from aiida import orm

In [2]: from ase.build import bulk
   ...: structure = orm.StructureData(ase=bulk('Co', 'hcp', 2.5, 4.06))
```

The `create_magnetic_configuration` can be used to quickly create a new `StructureData` with two different kinds:

```python
In [3]: from aiida_quantumespresso.calculations.functions.create_magnetic_configuration import create_magnetic_configuration
   ...:
   ...: results = create_magnetic_configuration(structure, [-2, 2])

In [4]: results['structure'].sites
Out[4]:
[<Site: kind name 'Co0' @ -2.7755575615629e-17,1.4433756729741,2.03>,
 <Site: kind name 'Co1' @ 0.0,0.0,0.0>]

In [5]: results['magnetic_moments'].get_dict()
Out[5]: {'Co0': 2, 'Co1': -2}
```

For more information, see the tutorial on how to work with magnetic systems:

https://aiida-quantumespresso.readthedocs.io/en/latest/tutorials/magnetism.html

The release also makes an important change in the dependencies related to a bug introduced in `pydantic`, see:

https://github.com/pydantic/pydantic/issues/5821

Hence the version of `pydantic` is adapted to `1.10.8`, where this bug is fixed.

### ‼️ Breaking changes

* `PwCalculation`: refactor `parent_folder` validation [[a389629](https://github.com/aiidateam/aiida-quantumespresso/commit/a389629387b74805ffe2f4d6515ac05b8f62b4d5)]

### ✨ New features

* Add the `create_magnetic_configuration` function [[d9b18a7](https://github.com/aiidateam/aiida-quantumespresso/commit/d9b18a7c20ce023018755c202f8d06cbf8bd27c5)]

### 👌 Improvements

* `PpParser`: Include exception in `ERROR_OUTPUT_DATAFILE_PARSE` message [[72f114e](https://github.com/aiidateam/aiida-quantumespresso/commit/72f114e4b05b45297abf0954b6334f5a461ed17e)]

### 🐛 Bug Fixes

* `PwParser`: Fix case when `settings` are not provided [[5d4a7d9](https://github.com/aiidateam/aiida-quantumespresso/commit/5d4a7d9405b757e2ecf65a72c1ee92aa2fb36a39)]

### 📚 Documentation

* Small improvements to "Installation" page [[90ad1d6](https://github.com/aiidateam/aiida-quantumespresso/commit/90ad1d6026d3c4b557970d6cc7626e85195ca4dc)]
* Switch to `sphinx-book-theme` from `pydata_sphinx_theme` [[3578a9d](https://github.com/aiidateam/aiida-quantumespresso/commit/3578a9d08af5e8cdd0851f852185bce2f2c7bd51)]
* Switch to using MyST Markdown [[37d2a14](https://github.com/aiidateam/aiida-quantumespresso/commit/37d2a14256480785491429c4ca424ebdc258ae34)]
* Add how-to for `PwCalculationTools.get_scf_accuracy` [[29b4db9](https://github.com/aiidateam/aiida-quantumespresso/commit/29b4db9e5c0d225331aba58981663ad4af641640)]
* Bump Python version for RTD build [[483d99b](https://github.com/aiidateam/aiida-quantumespresso/commit/483d99b77ca26ce4c607440922b432e39c16d1cf)]
* Fix breaking changes `CHANGELOG.md` [[7f4c4a1](https://github.com/aiidateam/aiida-quantumespresso/commit/7f4c4a108b09c31cfc9ee252d7e8fbe574ea477f)]

### 🔧 Maintenance

* Catch warning in `restart_mode` test for `PwBaseWorkChain` [[45e1907](https://github.com/aiidateam/aiida-quantumespresso/commit/45e19072f28d30fb4d6df1bbc489f38ce94cd1be)]
* Update `tox` configuration [[0702152](https://github.com/aiidateam/aiida-quantumespresso/commit/0702152c1dc18fd8d04cbf44f99a15761be955fe)]
* Add script to update `CHANGELOG.md` [[b21076c](https://github.com/aiidateam/aiida-quantumespresso/commit/b21076cdd1773fe3b7de18bc83e89ec7b367d837)]
* Remove Python 3.8 from the nightly workflow matrix [[0859d0a](https://github.com/aiidateam/aiida-quantumespresso/commit/0859d0a05c495c8a92208d118da6066f385600a7)]

### ⬆️ Update dependencies

* Restrict `pydantic` to at least `1.10.8` [[f430139](https://github.com/aiidateam/aiida-quantumespresso/commit/f430139e36723f73253efa66a6444018123760d2)]

### ♻️ Refactor

* `BasePwCpInputGenerator`: Remove superfluous validation [[0b6476c](https://github.com/aiidateam/aiida-quantumespresso/commit/0b6476c6c8379a0ea6575f4323979d136f7220aa)]
*  Move basic parsing into `BaseParser` class [[1c223c7](https://github.com/aiidateam/aiida-quantumespresso/commit/1c223c78f0a2eda38089a08d9f21626f49d2fd6b)]

## v4.3.0

Release version `4.3.0` comes with a lot of new features, improvements and bug fixes.
Although this is technically a minor release, this new version does come with some _minor_ breaking changes:

* The default value of `clean_workdir` was changed from `True` to `False`.
* The `automatic_parallelization` feature of the `PwBaseWorkChain` was removed, as this was badly broken with no clear path to fixing it. Moreover, [v7.1 of Quantum ESPRESSO](https://gitlab.com/QEF/q-e/-/releases/qe-7.1) has implemented some basic automated parallelization for `pw.x` when no parallelization flags are specified.
* Work chains no longer set/override static inputs inside the steps of their outline. All defaults are moved to the protocol files instead, and it is recommended to always use the `get_builder_from_protocol()` method to run a work chain. For more details, [see the corresponding PR](https://github.com/aiidateam/aiida-quantumespresso/pull/902).

Support for Quantum ESPRESSO v7.2 `pw.x` parsing was added, as well as the `HubbardStructureData` data plugin and corresponding implementation in the `PwCalculation` plugin to support the [new input syntax for using Hubbard corrections](https://gitlab.com/QEF/q-e/-/releases/qe-7.1#incompatible-changes-in-71-version) in the `pw.x` code.

For `ph.x`, the symmetry labels printed in the `stdout` after the mode frequencies are now parsed for each q-point.
In the case of restarts in the `PhBaseWorkChain`, care has been taken to properly parse these symmetry labels from the separate output files and merge them, so the `output_parameters` are the same with or without restarts.

The `XpsWorkChain` has been added to compute the XPS spectra using core-hole pseudo potentials and the `pw.x` code. The work chain can handle both molecules and crystals. The chemical shifts, as well as the cole-level binding energy can be obtained.

Finally, improved error handling for the `PwBaseWorkChain` was implemented for several typical failure modes of the `pw.x` code, and the `CRASH` file is now also retrieved, as this will be the default location of error messages from Quantum ESPRESSO v7.2 onwards.

### ‼️ Breaking changes

* Protocols: Set `clean_workdir` default to `False` [[f8e512f](https://github.com/aiidateam/aiida-quantumespresso/commit/f8e512f9cb5404d702aa1d721e2369c7257f977f)]
* `PwBaseWorkChain`: Remove `automatic_parallelization` input [[5cae75f](https://github.com/aiidateam/aiida-quantumespresso/commit/5cae75ff671268dc1e5684e1c84f801a10839e0b)]
* Protocols: Move all static work chain inputs to protocol [[01f1470](https://github.com/aiidateam/aiida-quantumespresso/commit/01f14701e6846338f84db25bf7e654fcdfb6f927)]

### ✨ New features

* `PwBaseWorkChain`: Add `ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY` handler [[9291f84](https://github.com/aiidateam/aiida-quantumespresso/commit/9291f841445800fa2b38d97b6957e98a2818b624)]
* `PwCalculation`: Add exit code and handler for `scale_h` error [[d350d7e](https://github.com/aiidateam/aiida-quantumespresso/commit/d350d7eb5b7bfcf4682ee8f635c721708be1d1d7)]
* `PwParser`: Add retrieval and parsing of `CRASH` file [[7f53c96](https://github.com/aiidateam/aiida-quantumespresso/commit/7f53c96bf744ed622fe7b1d13b0f7e9198ed5656)]
* Add the `HubbardStructureData` data plugin [[355020c](https://github.com/aiidateam/aiida-quantumespresso/commit/355020ca9ea63db0803ecc5746f93a879f488696)]
* `PwParser`: Add support for QE v7.2 [[57f5f8f](https://github.com/aiidateam/aiida-quantumespresso/commit/57f5f8ffa359ed5f5d937c43ae1f53cbc55eb314)]
* `PhParser`: parse symmetry labels from the `stdout` of `ph.x` [[8a9950a](https://github.com/aiidateam/aiida-quantumespresso/commit/8a9950aa870a26badb24b572910f7b8526236f37)]
* Add Feature: `XpsWorkChain` [[a9d124e](https://github.com/aiidateam/aiida-quantumespresso/commit/a9d124e2cfa5a5a5891affbe40515d38d1dd3913)]
* Add Feature: `XspectraCrystalWorkChain` [[01e7593](https://github.com/aiidateam/aiida-quantumespresso/commit/01e759319552645072b7a89dd0347225ab14d635)]
* Add `OpenGridCalculation` [[361ff04](https://github.com/aiidateam/aiida-quantumespresso/commit/361ff0413659b23f895747e04f8448f574811958)]
* Feature/get xspectra structures [[bc63d56](https://github.com/aiidateam/aiida-quantumespresso/commit/bc63d5661a6e79e3cf599ada3fdfc8eb3cb8193d)]

### 👌 Improvements

* ``PwBaseWorkChain``: improve diagonalization handler [[cc29488](https://github.com/aiidateam/aiida-quantumespresso/commit/cc29488c9e450bac5af7d28ee3f79ab8bcb01b6a)]
* `PwCalculation`: Add `ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY` code [[92a6c6f](https://github.com/aiidateam/aiida-quantumespresso/commit/92a6c6fb00c02ace0f8afc0d06d8f938a0e2d221)]
* `PwParser`: Keep scheduler parser error unless more specific available [[b8d6a3a](https://github.com/aiidateam/aiida-quantumespresso/commit/b8d6a3a489c8d24a06cfc58956a97b121ad92316)]
* Protocols: consider `pbc` of `StructureData` [[ef21642](https://github.com/aiidateam/aiida-quantumespresso/commit/ef21642ba8c715b384e726628f866f7a3df63890)]
* CLI: Simplify `cli.utils.validate_hubbard_parameters` and add tests [[74d25d1](https://github.com/aiidateam/aiida-quantumespresso/commit/74d25d1782cddb2801781999e6c708191f1cdae7)]
* `PwCalculation`: Add exit codes, mostly related to diagonalization [[159b83e](https://github.com/aiidateam/aiida-quantumespresso/commit/159b83e5a96db258697fa821b48746adfbacd2ad)]
* `PwParser`: Include path in exception when XML schema not found [[56fdb57](https://github.com/aiidateam/aiida-quantumespresso/commit/56fdb57643fa44e7bcff82640d9cbcb26f7759b3)]

### 🐛 Bug Fixes

* `PwCalculation`: Fix restart validation for `nscf`/`bands` [[29a0dfa](https://github.com/aiidateam/aiida-quantumespresso/commit/29a0dfa3257863f73de17a0347dfa578da5d292e)]
* `PwParser`: Correct the `disk_io` type for all XML schemas [[89d39a4](https://github.com/aiidateam/aiida-quantumespresso/commit/89d39a47d4eae0f9b39ef93fe07404a5102c6ee7)]
* `PwParser`: Fix trajectory verification for `fixed_coords` [[105670f](https://github.com/aiidateam/aiida-quantumespresso/commit/105670fd145c572213bb750069679ef6fca7a80d)]
* `XspectraCoreWorkChain`: Ensure `overrides` respected in `get_builder_from_protocol()` [[7b2d701](https://github.com/aiidateam/aiida-quantumespresso/commit/7b2d701cb7f61719efa5a14f762b8f85078fb4dd)]
* `PhCalculation`: Fix exception when `parent_folder` is on different computer [[24b3779](https://github.com/aiidateam/aiida-quantumespresso/commit/24b3779e0df2d30378898f89dd884b69593218e9)]

### 📚 Documentation

* Fix the intersphinx URL of the AiiDA documentation [[08532d8](https://github.com/aiidateam/aiida-quantumespresso/commit/08532d8f69d87bd1853bb2a9010854efbc6bbaee)]
* `README.md`: Fix Build Status badge to point to `ci` action [[1e9c345](https://github.com/aiidateam/aiida-quantumespresso/commit/1e9c3454174807e115f48eb7667f928680a11dfc)]
* Add developer section to README [[146bc57](https://github.com/aiidateam/aiida-quantumespresso/commit/146bc57d74b76a11d7205a0f8449ae39c7675110)]

### 🔧 Maintenance

* DevOps: Symlink `CRASH` and `data-file.xml` files, and reduce fixture file contents [[0a48533](https://github.com/aiidateam/aiida-quantumespresso/commit/0a48533f40e276e5423428e0ca9745bd0b789121)]

### ⬆️ Update dependencies
* Symlink `CRASH` and `data-file.xml` files, and reduce fixture file contents [[0a48533](https://github.com/aiidateam/aiida-quantumespresso/commit/0a48533f40e276e5423428e0ca9745bd0b789121)]
* Devops: Update `pylint` version [[9f4142d](https://github.com/aiidateam/aiida-quantumespresso/commit/9f4142d4713b0d1c32042b3ca35905d725954bf4)]

### ♻️ Refactor
* Update `pylint` version [[9f4142d](https://github.com/aiidateam/aiida-quantumespresso/commit/9f4142d4713b0d1c32042b3ca35905d725954bf4)]


## v4.2.0

### Features
- Add the `XspectraBaseWorkChain` [[#854]](https://github.com/aiidateam/aiida-quantumespresso/pull/854)
- Add the `XspectraCoreWorkChain` [[#872]](https://github.com/aiidateam/aiida-quantumespresso/pull/872)

### Changes
- Protocols: Use v1.2 of SSSP and PBEsol by default for pw.x protocols [[#884]](https://github.com/aiidateam/aiida-quantumespresso/pull/884)

### Fixes
- `PwBandsWorkChain`: Fix `nbnd` for spin-polarised calculations [[#875]](https://github.com/aiidateam/aiida-quantumespresso/pull/875)

### Dependencies
- Update requirement `aiida-pseudo~=1.0` [[#882]](https://github.com/aiidateam/aiida-quantumespresso/pull/882)
- Update pre-commit requirement `isort==5.12.0` [[#883]](https://github.com/aiidateam/aiida-quantumespresso/pull/883)

### Devops
- Fix nightly build tests [[#870]](https://github.com/aiidateam/aiida-quantumespresso/pull/870)

### Documentation
- Add various `settings` parameters for `PwCalculation` [[#869]](https://github.com/aiidateam/aiida-quantumespresso/pull/869)


## v4.1.0

### New

- Add the `options` argument to all `get_builder_from_protocol` [[#846]](https://github.com/aiidateam/aiida-quantumespresso/pull/846)
- `XspectraParser`: Parse parameters from `CalcJob` inputs [[#853]](https://github.com/aiidateam/aiida-quantumespresso/pull/853)

### Fixes

- `XspectraCalculation`/`XspectraParser`: various improvements and fixes [[#837]](https://github.com/aiidateam/aiida-quantumespresso/pull/837)
- CLI: update to be compatible with `aiida-core==2.1` [[#855]](https://github.com/aiidateam/aiida-quantumespresso/pull/855)
- Fix bug in `wrap_bare_dict_inputs` [[#860]](https://github.com/aiidateam/aiida-quantumespresso/pull/860)
- `PwBaseWorkChain`: Fix kpoints override in `get_builder_from_protocol` [[#862]](https://github.com/aiidateam/aiida-quantumespresso/pull/862)
- `BasePwCpInputGenerator`: do not require pseudos in `validate_inputs` [[#866]](https://github.com/aiidateam/aiida-quantumespresso/pull/866)
- `PwBaseWorkChain`: Set `restart_mode` in `parameters` if `parent_folder` [[#867]](https://github.com/aiidateam/aiida-quantumespresso/pull/867)

### Deprecation

- Remove obsolete and broken TCOD related code [[#850]](https://github.com/aiidateam/aiida-quantumespresso/pull/850)

### Devops

- DevOps: Add a job with integration tests to the CI workflow

### Dependencies

- Dependencies: Add support for Python 3.11 [[#857]](https://github.com/aiidateam/aiida-quantumespresso/pull/857)


## v4.0.1

### Fixes
- `PwCalculation/PhCalculation`: Add validation for `pseudos` port [[#839]](https://github.com/aiidateam/aiida-quantumespresso/pull/839)
- `Q2rParser`: Set `filename` for `force_constants` [[#827]](https://github.com/aiidateam/aiida-quantumespresso/pull/827)
- Docs: Update the `Get Started` example for `PhCalculation` [[#829]](https://github.com/aiidateam/aiida-quantumespresso/pull/829)
- DevOps: Add `AIIDA_WARN_V3` environment variable to CI and address all deprecation warnings [[#836]](https://github.com/aiidateam/aiida-quantumespresso/pull/836)

### Dependencies
- Add lower requirement `aiida-core>=2.0.4`, fixes error in `NamelistCalculation` restarting from `FolderData` [[#841]](https://github.com/aiidateam/aiida-quantumespresso/pull/841)


## v4.0.0

This major release mainly provides support for `aiida-core` v2.0, as well as the `xspectra.x` code of Quantum ESPRESSO.
Version 7.1 of Quantum ESPRESSO is also now supported.
Note that in line with our new [compatibility policy](https://github.com/aiidateam/aiida-quantumespresso#compatibility-matrix), from this major version onwards we will no longer be supporting Quantum ESPRESSO versions below v6.6.

### New

- `PhBaseWorkChain`: Add a first draft for the protocols [[#810]](https://github.com/aiidateam/aiida-quantumespresso/pull/810)
- Add the `XspectraCalculation` and `XspectraParser` plugins [[#815]](https://github.com/aiidateam/aiida-quantumespresso/pull/815)

### Improvements

- `PwCalculation`: allow regex patterns in detect_important_message [[#771]](https://github.com/aiidateam/aiida-quantumespresso/pull/771)
- Parsers: add exception title to `ERROR_UNEXPECTED_PARSER_EXCEPTION` [[#793]](https://github.com/aiidateam/aiida-quantumespresso/pull/793)
- `CalcJob`: simplify `parent_folder` handling [[#805]](https://github.com/aiidateam/aiida-quantumespresso/pull/805)
- `MatDynCalculation`: use validator to check `flfrc` input [[#814]](https://github.com/aiidateam/aiida-quantumespresso/pull/814)
- `CpParser`: add stresses and fix for QE v7.0 [[#786]](https://github.com/aiidateam/aiida-quantumespresso/pull/786)
- `PwParser`: Add support for QE v7.1 [[#822]](https://github.com/aiidateam/aiida-quantumespresso/pull/822)

### Bug Fixes

- `PhCalculation`: Fix bug when `only_initialization=True` [[#813]](https://github.com/aiidateam/aiida-quantumespresso/pull/813)

### Devops

- Devops: fix the nightly build [[#784]](https://github.com/aiidateam/aiida-quantumespresso/pull/784)
- DevOps: Move the source directory into `src/` [[#812]](https://github.com/aiidateam/aiida-quantumespresso/pull/812)
- DevOps: update nightly to show as failed on error [[#817]](https://github.com/aiidateam/aiida-quantumespresso/pull/817)
- DevOps: only run Slack notification if tests failed [[#819]](https://github.com/aiidateam/aiida-quantumespresso/pull/819)
- DevOps: address deprecation warnings [[#820]](https://github.com/aiidateam/aiida-quantumespresso/pull/820)

### Dependencies

- Dependencies: put upper limit markupsafe<2.1 [[#790]](https://github.com/aiidateam/aiida-quantumespresso/pull/790)
- Dependencies: Update Python requirements [[#803]](https://github.com/aiidateam/aiida-quantumespresso/pull/803)
- Dependencies: Update requirement for `aiida-core~=2.0` [[#804]](https://github.com/aiidateam/aiida-quantumespresso/pull/804)
- Refactor: Additional compatibility for `aiida-core` v2.0 [[#806]](https://github.com/aiidateam/aiida-quantumespresso/pull/806)
- Remove the use of deprecated distutils package [[#808]](https://github.com/aiidateam/aiida-quantumespresso/pull/808)


## v3.5.2

### New

- `PwParser`: Add support for QE v7.0 [[#774]](https://github.com/aiidateam/aiida-quantumespresso/pull/774)

### Bug Fixes

- Fix XML parsing for structures with numbered kinds [[#770]](https://github.com/aiidateam/aiida-quantumespresso/pull/770)

### Documentation

- Update pw_tutorial.rst [[#769]](https://github.com/aiidateam/aiida-quantumespresso/pull/769)

### Devops

- Add GitHub Actions workflow for continuous deployment [[#776]](https://github.com/aiidateam/aiida-quantumespresso/pull/776)
- Adopt PEP 621 and move build spec to pyproject.toml [[#779]](https://github.com/aiidateam/aiida-quantumespresso/pull/779)
- CI: update the GHA ci.yml workflow [[#780]](https://github.com/aiidateam/aiida-quantumespresso/pull/780)
- Update pre-commit dependencies and configuration [[#781]](https://github.com/aiidateam/aiida-quantumespresso/pull/781)
- CI: add the isort pre-commit hook [[#781]](https://github.com/aiidateam/aiida-quantumespresso/pull/781)


## v3.5.1

### Improvements
- `PhBaseWorkChain`: add handler for diagonalization errors [[#757]](https://github.com/aiidateam/aiida-quantumespresso/pull/757)
- `PhBaseWorkChain`: add handler for `ERROR_SCHEDULER_OUT_OF_WALLTIME` [[#754]](https://github.com/aiidateam/aiida-quantumespresso/pull/754)

### Documentation
- Fix "suggest edit" link [[#768]](https://github.com/aiidateam/aiida-quantumespresso/pull/768)

### Packaging
- Add `calculations/helpers/*.xml` to `MANIFEST.in` [[#760]](https://github.com/aiidateam/aiida-quantumespresso/pull/760)


## v3.5.0

In this minor version release we provide support for Quantum ESPRESSO v6.8 and add several improvements and bug fixes mostly related to input overrides and restarting calculations.
Another important addition is the introduction of [a compatibility policy for Quantum ESPRESSO](https://github.com/aiidateam/aiida-quantumespresso/blob/develop/README.md#compatibility-matrix), which will go into effect starting v4.0.0.

### New

- `PwParser`: add support for QE v6.8 with new XML schema file [[#717]](https://github.com/aiidateam/aiida-quantumespresso/pull/717)
- `PwBaseWorkChain`: handle diagonalization issues [[#744]](https://github.com/aiidateam/aiida-quantumespresso/pull/744)

### Improvements

- `PwBaseWorkChain`: add handler for intentional non-convergence [[#695]](https://github.com/aiidateam/aiida-quantumespresso/pull/695)
- `PwBaseWorkChain`: make magnetism from `overrides` absolute [[#731]](https://github.com/aiidateam/aiida-quantumespresso/pull/731)
- `PwBaseWorkChain`: Improve restart and validate inputs [[#722]](https://github.com/aiidateam/aiida-quantumespresso/pull/722)
- `PwRelaxWorkChain`: do not override SYSTEM/nbnd in final scf if set. [[#708]](https://github.com/aiidateam/aiida-quantumespresso/pull/708)
- `PwCalculation`: adjust parent_folder validation [[#741]](https://github.com/aiidateam/aiida-quantumespresso/pull/741)
- `PwBaseWorkChain`: make overrides absolute [[#741]](https://github.com/aiidateam/aiida-quantumespresso/pull/741)
- `projwfc.x`: parse from XML instead of parent calc [[#747]](https://github.com/aiidateam/aiida-quantumespresso/pull/747)

### Bug Fixes

- `PhCalculation`: use `verbosity` instead of `iverbosity` [[#633]](https://github.com/aiidateam/aiida-quantumespresso/pull/633)
- Protocols: Add settings from `overrides` [[#725]](https://github.com/aiidateam/aiida-quantumespresso/pull/725)
- `PhBaseWorkChain`: spec exposed outputs of incorrect process class [[#732]](https://github.com/aiidateam/aiida-quantumespresso/pull/732)

### Dependencies

- Temporarily set upper limit for `psycopg2` [[#707]](https://github.com/aiidateam/aiida-quantumespresso/pull/707)

### Documentation

- Docs: add compatibility policy for Quantum ESPRESSO [[#737]](https://github.com/aiidateam/aiida-quantumespresso/pull/737)
- `README.md`: cleanup the syntax of compatibility matrix [[#738]](https://github.com/aiidateam/aiida-quantumespresso/pull/738)
- `README.md`: add shields for Quantum ESPRESSO compatibility [[#739]](https://github.com/aiidateam/aiida-quantumespresso/pull/739)


## v3.4.2

This patch release fixes several bugs in the protocols feature, adds support for the `CutoffsPseudoPotentialFamily` introduced in `aiida-pseudo==0.6.1` and adapts the `ProtocolMixin` class so it can also be used by other packages.
The legacy `PwBandStructureWorkChain` is now also fully deprecated and will be removed in the next major release.
It has been replaced by the `PwBandsWorkChain`, which is now equipped with its own protocol through the `get_builder_from_protocol()` method.

### Improvements
- `ProtocolMixin`: increase flexibility and usability for other packages [[#678]](https://github.com/aiidateam/aiida-quantumespresso/pull/678)
- Protocols: Add support for `CutoffsPseudoPotentialFamily` [[#684]](https://github.com/aiidateam/aiida-quantumespresso/pull/684)

### Bug Fixes
- Protocols: Add usage of metadata and parallelisation `overrides` [[#652]](https://github.com/aiidateam/aiida-quantumespresso/pull/652)
- `PwBaseWorkChain`: fix bug in `validate_resources` validator [[#683]](https://github.com/aiidateam/aiida-quantumespresso/pull/683)

### Deprecation
- `PwBandStructureWorkChain`: Deprecate and fix [[#688]](https://github.com/aiidateam/aiida-quantumespresso/pull/688)


## v3.4.1

### Dependencies
- Dependencies: update to `aiida-pseudo==0.6.0` [[#671]](https://github.com/aiidateam/aiida-quantumespresso/pull/#671)

### Documentation
- Docs: fix `gamma_only` value in `settings_dict` [[#672]](https://github.com/aiidateam/aiida-quantumespresso/pull/#672)

### Devops
- CI: use `postgres-12` for Github Actions since `ubuntu-latest` was updated to Ubuntu Focal 20.04 [[#674]](https://github.com/aiidateam/aiida-quantumespresso/pull/#674)


## v3.4.0

### Features
Add the `PdosWorkChain` (#418)
Add support for `pw.x` and `cp.x` v6.7 (#626)
`PwCalculation`: Add parsing of up/dw Fermi energy (#622)

### Changes
Add support for Python 3.9 (#666)
`PwRelaxWorkChain`: remove compatibility for volume optimization (only affects the input generation protocols added in beta phase in `v3.3.0`) (#660)
`RelaxType`: Rename `ATOMS` to `POSITIONS`(only affects the input generation protocols added in beta phase in `v3.3.0`)  (#658)


## v3.3.1

### Bug fixes
- Protocols: fix unintended mutation of `overrides` dictionary passed as argument [[#650]](https://github.com/aiidateam/aiida-quantumespresso/pull/#650)
- Protocols: fix ionic convergence thresholds for `fast` protocol [[#642]](https://github.com/aiidateam/aiida-quantumespresso/pull/#642)

### Dependencies
- Update requirement to `aiida-pseudo~=0.5.0` [[#646]](https://github.com/aiidateam/aiida-quantumespresso/pull/#646)

### Devops
- Fix nightly build vs `aiida-core` develop branch by dropping Python 3.6 and install `aiida-core` before `aiida-quantumespresso` [[#646]](https://github.com/aiidateam/aiida-quantumespresso/pull/#646)
- Fix documentation and add tox configuration [[#643]](https://github.com/aiidateam/aiida-quantumespresso/pull/#643)


## v3.3.0
The input generation protocols are a new feature that are now in a beta stage.
The API might still change in future releases, and it has not been fully tested.

### Features
- **BETA** Add input generation protocols for the pw.x work chains [[#608]](https://github.com/aiidateam/aiida-quantumespresso/pull/608), [[#616]](https://github.com/aiidateam/aiida-quantumespresso/pull/616), [[#631]](https://github.com/aiidateam/aiida-quantumespresso/pull/631), [[#632]](https://github.com/aiidateam/aiida-quantumespresso/pull/631), [[#634]](https://github.com/aiidateam/aiida-quantumespresso/pull/634), [[#635]](https://github.com/aiidateam/aiida-quantumespresso/pull/635), [[#638]](https://github.com/aiidateam/aiida-quantumespresso/pull/638)
- `CpCalculation`: add support for `AUTOPILOT` mode [[#455]](https://github.com/aiidateam/aiida-quantumespresso/pull/455)
- `PwCalculation`: add explicit `parallelization` input port [[#554]](https://github.com/aiidateam/aiida-quantumespresso/pull/554)
- `ProjwfcParser`: support kind names with numerals, dashes and underscores [[#591]](https://github.com/aiidateam/aiida-quantumespresso/pull/591)

### Improvements
- `PwRelaxWorkChain`: expose `PwBaseWorkChain` separately for final SCF [[#569]](https://github.com/aiidateam/aiida-quantumespresso/pull/569)
- Replace old format string interpolation with f-strings [[#579]](https://github.com/aiidateam/aiida-quantumespresso/pull/579)
- `Parsers`: remove explicit check for the retrieved folder's existence [[#580]](https://github.com/aiidateam/aiida-quantumespresso/pull/580)
- `BasePwCpInputGenerator`: remove duplicate setting cmdline arguments [[#585]](https://github.com/aiidateam/aiida-quantumespresso/pull/585)
- `PpParser`: reduce memory usage by deleting raw data once parsed [[#587]](https://github.com/aiidateam/aiida-quantumespresso/pull/587)
- `ForceConstantsData`: add `filename` and `**kwargs` to constructor [[#599]](https://github.com/aiidateam/aiida-quantumespresso/pull/599)
- `PwBaseWorkChain`: restart from scratch when convergence not reached [[#604]](https://github.com/aiidateam/aiida-quantumespresso/pull/604)
- `PwParser`: ignore stress threshold if `CELL.cell_dofree != 'all'` [[#613]](https://github.com/aiidateam/aiida-quantumespresso/pull/613)
- `PwBaseWorkChain`: add report on disabling bands sanity check [[#629]](https://github.com/aiidateam/aiida-quantumespresso/pull/629)

### Deprecated
- `PwBaseWorkChain`: deprecate the `pseudo_family` input [[#603]](https://github.com/aiidateam/aiida-quantumespresso/pull/603)
- `PwRelaxWorkChain`: replace `relaxation_scheme` with `relax_type` [[#614]](https://github.com/aiidateam/aiida-quantumespresso/pull/614)
- `PwRelaxWorkChain`: deprecate the `final_scf` input [[#569]](https://github.com/aiidateam/aiida-quantumespresso/pull/569)

### Removed
- Remove parser information from `output_parameters` [[#597]](https://github.com/aiidateam/aiida-quantumespresso/pull/597)

### Bug Fixes
- Fix minor bug in `utils.protocols._get_all_protocol_modifiers` [[#592]](https://github.com/aiidateam/aiida-quantumespresso/pull/592)
- `PwCalculation`: fix bug when `hubbard_file` is specified [[#596]](https://github.com/aiidateam/aiida-quantumespresso/pull/596)
- `ProjwfcParser`: fix bug in case of for more than 1000 orbitals [[#624]](https://github.com/aiidateam/aiida-quantumespresso/pull/624)

### Dependencies
- Drop support for Python 3.5 [[#578]](https://github.com/aiidateam/aiida-quantumespresso/pull/578)
- Add support for `aiida-pseudo` [[#581]](https://github.com/aiidateam/aiida-quantumespresso/pull/581), [[#630]](https://github.com/aiidateam/aiida-quantumespresso/pull/630)
- Update the `aiida-pseudo` requirements [[#609]](https://github.com/aiidateam/aiida-quantumespresso/pull/609)

### Devops
- Add the `Framework :: AiiDA` trove classifier to the package [[#584]](https://github.com/aiidateam/aiida-quantumespresso/pull/584)
- CLI: add tests for the process launch commands [[#589]](https://github.com/aiidateam/aiida-quantumespresso/pull/589)
- Pre-commit: re-enable formatting after `define` methods [[#619]](https://github.com/aiidateam/aiida-quantumespresso/pull/619)


## v3.2.1:

### Improvements
- `PpParser`: reduce memory usage by deleting raw data once parsed [[#587]](https://github.com/aiidateam/aiida-quantumespresso/pull/587)
- `ProjwfcParser`: support kind names with numerals, dashes and underscores [[#591]](https://github.com/aiidateam/aiida-quantumespresso/pull/591)


## v3.2.0:

### Features
- `Pw/CpCalculation`: allow explicitly specifying `ibrav` [[#503]](https://github.com/aiidateam/aiida-quantumespresso/pull/503))
- `PwParser`: add support for output format of pw.x v6.6 [[#552]](https://github.com/aiidateam/aiida-quantumespresso/pull/552)
- `ProjwfcParser`: add support for non-collinear and spinorbit calculations [[#546]](https://github.com/aiidateam/aiida-quantumespresso/pull/546)
- `NamelistCalculation`: allow symlinking `parent_folder` instead of copying [[#555]](https://github.com/aiidateam/aiida-quantumespresso/pull/555)

### Bug fixes
- `PwParser`: fix bug in raw parsing of van der Waals and stress [[#550]](https://github.com/aiidateam/aiida-quantumespresso/pull/550)
- `PwBaseWorkChain`: do no assume `output_band` exists in band sanity check [[#544]](https://github.com/aiidateam/aiida-quantumespresso/pull/544)
- `PwBandsWorkChain`: properly mark optional outputs as not required [[#571]](https://github.com/aiidateam/aiida-quantumespresso/pull/571)
- `PwBandsWorkChain`: incorrect number of bands for non-collinear spin [[#570]](https://github.com/aiidateam/aiida-quantumespresso/pull570/)
- `PwBandStructureWorkChain`: add context argument to validator [[#560]](https://github.com/aiidateam/aiida-quantumespresso/pull/560)
- Fix the incorrect usage of repository interface of `aiida-core` [[#549]](https://github.com/aiidateam/aiida-quantumespresso/pull/549))

### Dependencies
- Update to qe-tools 2.0.0rc1. [[#529]](https://github.com/aiidateam/aiida-quantumespresso/pull/529)
- Increase minimum required `aiida-core` version to 1.3 [[#573]](https://github.com/aiidateam/aiida-quantumespresso/pull/573))

### Devops
- CI: add nightly build against `aiida-core/develop` [[#568]](https://github.com/aiidateam/aiida-quantumespresso/pull/568)
- Pre-commit: move `pytest` and `pylint` config to `pyproject.toml` [[#562]](https://github.com/aiidateam/aiida-quantumespresso/pull/562)


## v3.1.0:

### Changes
- Drop Python 2 support [[#515]](https://github.com/aiidateam/aiida-quantumespresso/pull/515)
- Move calcfunctions to `aiida_quantumespresso.calculations.functions` [[#520]](https://github.com/aiidateam/aiida-quantumespresso/pull/520)
- Replace builtin `BaseRestartWorkChain` with the one from `aiida-core` [[#519]](https://github.com/aiidateam/aiida-quantumespresso/pull/519)

### Features
- `PpCalculation`: reinstate the `settings` input node [[#537]](https://github.com/aiidateam/aiida-quantumespresso/pull/537)
- `PpCalculation`: add support for retrieving and parsing multiple files [[#533]](https://github.com/aiidateam/aiida-quantumespresso/pull/533)

### Bug fixes
- `PpParser`: improve the Gaussian cube file parsing algorithm. [[#535]](https://github.com/aiidateam/aiida-quantumespresso/pull/535)
- `PpParser`: fix whitespace bug in the `parse_gnuplot2D` method [[#534]](https://github.com/aiidateam/aiida-quantumespresso/pull/534)
- Fix import or `URLError` which has been removed in `xmlschema` [[#524]](https://github.com/aiidateam/aiida-quantumespresso/pull/524)


## v3.0.0:
This is the first version to support `aiida-core~=1.0` and since the API of `aiida-core` changed significantly, the code of `aiida-quantumespresso` did as well.
This means that `aiida-quantumespresso` `v2.0` and `v3.0` are fundamentally incompatible, which is why we do not include an explicit change log list.


## v2.1.0:
Minor release adding a new feature, fixing a critical bug in `PwParser` and reducing the size of the package

*Nota Bene*: due to a critical bug discoverd in the `PwParser` (see below) the parser version has been updated to `v2.1.0`.
Calculations parsed with an older version potentially contain incorrect reduced symmetries in the `output_parameters` node.
Instructions to check the version of the parser that was used, [can be found in the documentation](http://aiida-quantumespresso.readthedocs.io/en/latest/user_guide/calculation_plugins/pw.html#parser-version)

### Improvements
- Added support for the `ATOMIC_FORCES` card in the `PwCalculation` plugin [[#168]](https://github.com/aiidateam/aiida-quantumespresso/pull/168)
- Tests and fixture data have been removed from the main package, massively reducing distribution size [[#154]](https://github.com/aiidateam/aiida-quantumespresso/pull/154)

### Critical bug fixes
- The mapping of the raw symmetry operations onto a reduced set in the `PwParser` contained a bug and mapped incorrect rotations [[#171]](https://github.com/aiidateam/aiida-quantumespresso/pull/171)

### Minor bug fixes
- Fix a bug in `ProjwfcParser` where due to improper sorting, projections were associated with the wrong output files [[#165]](https://github.com/aiidateam/aiida-quantumespresso/pull/165)
- Adapt the keyword `type` to `node_type` in the `Node.get_outputs` method [[#143]](https://github.com/aiidateam/aiida-quantumespresso/pull/143)
- Fix a missing parameter in the string formatter of one of the error handlers of `PwBaseWorkChain` [[#158]](https://github.com/aiidateam/aiida-quantumespresso/pull/158)


## v2.0.1:
Patch release with some small bug fixes

### Update to newer versions of Quantum ESPRESSO
- Add the input helper XML files for v6.0, v6.1 and v6.2 [[#135]](https://github.com/aiidateam/aiida-quantumespresso/pull/135)

### Bug fixes
- Parse the standard output in the `PwParser` even if XML file is missing in order to still get errors in output parameters [[#133]](https://github.com/aiidateam/aiida-quantumespresso/pull/133)
- Bugfix for the `_error_handlers` list attribute from `BaseRestartWorkChain` that is now appended to the correct list[[#127]](https://github.com/aiidateam/aiida-quantumespresso/pull/127)
- Number of minor bugfixes to the `PwParser` to adapt to the updated AiiDA API [[#123]](https://github.com/aiidateam/aiida-quantumespresso/pull/123) [[#138]](https://github.com/aiidateam/aiida-quantumespresso/pull/138)
- Fix a bug causing cmdline args to be ignored in `namelists.py` [[#122]](https://github.com/aiidateam/aiida-quantumespresso/pull/122)


## v2.0.0:
Major release with a lot of new functionality, mostly centered around the workflows

### Improvements
- Implemented the `BaseRestartWorkChain` which defines much of the required scaffolding for any base workchain that launches a calculation [[#75]](https://github.com/aiidateam/aiida_core/pull/75)
- Add the `Q2rBaseWorkChain` and `MatdynBaseWorkChain` [[#93]](https://github.com/aiidateam/aiida_core/pull/93)
- Use the update wrappers for `SeeKpath` in `aiida-core` in the `PwBandsWorkChain` [[#104]](https://github.com/aiidateam/aiida_core/pull/104)
- Removed the use of `deepcopy` in all workchains, which can lead to bugs or unwanted duplication of nodes [[#118]](https://github.com/aiidateam/aiida_core/pull/118)

### Backwards incompatible changes
- Command line interface has been migrated to use `click` and are located in `aiida.cli` [[#81]](https://github.com/aiidateam/aiida_core/pull/81)
- `PhBaseWorkChain`: input `parent_calc` was renamed to `parent_folder` [[#87]](https://github.com/aiidateam/aiida_core/pull/87)

### Calculations
- `PwCalculation`: bands are now parsed by default making use of the new `retrieve_temporary_file_list` concept in `aiida-core` [[#36]](https://github.com/aiidateam/aiida_core/pull/36)
- `PwCalculation`: new option for `settings` to parse the atomic occupations [[#55]](https://github.com/aiidateam/aiida_core/pull/55)
- `PwCalculation`: add parsing of the electronic and ionic dipole if `lelfield` is used [[#105]](https://github.com/aiidateam/aiida_core/pull/105)

### Parsers
- `Pw2wannier90Parser`: added the parser for the `Pw2wannier90Calculation` class [[#38]](https://github.com/aiidateam/aiida_core/pull/38)

### Workflows
- `PwBaseWorkChain`: added `output_band` as optional output node [[#29]](https://github.com/aiidateam/aiida_core/pull/29)
- `PwBaseWorkChain`: added `output_array` as optional output node [[#97]](https://github.com/aiidateam/aiida_core/pull/97)
- `PwBaseWorkChain`: implemented automatic parallelization [[#39]](https://github.com/aiidateam/aiida_core/pull/39)
- `PhBaseWorkChain`: new optional input `only_initialization` to run an initialization calculation [[#101]](https://github.com/aiidateam/aiida_core/pull/101)
- `PwRelaxWorkChain`: new optional input `kpoints_force_parity` for generation of kpoint mesh [[#61]](https://github.com/aiidateam/aiida_core/pull/61)
- `PwRelaxWorkChain`: new optional input `max_meta_convergence_iterations` to limit the number of volume convergence steps [[#66]](https://github.com/aiidateam/aiida_core/pull/66)
- `PwBandsWorkChain`: make relaxing of structure optional [[#46]](https://github.com/aiidateam/aiida_core/pull/46)

### Bug fixes
- `PwRelaxWorkChain`: fixed bug when both `clean_workdir` and `final_scf` were enabled [[#59]](https://github.com/aiidateam/aiida_core/pull/59)
- `PwRelaxWorkChain`: properly unwrap the value of the `relaxation_scheme` input in the parameter input node of `PwCalculation` [[#119]](https://github.com/aiidateam/aiida_core/pull/119)
- `PwBaseWorkChain`: correctly set `do_break` for the `ErrorHandlerReport` of certain error handlers [[#113]](https://github.com/aiidateam/aiida_core/pull/113)


## v1.0.1:
Minor patch release with some small bug fixes

### Minor bug fixes
- Fix entry point for the `ForceconstantsData` class [[#22]](https://github.com/aiidateam/aiida_core/pull/22)
- Fix entry point for the `PwimmigrantCalculation` class [[#24]](https://github.com/aiidateam/aiida_core/pull/24)


## v1.0.0:
First official release of `aiida-quantumespresso`, the official plugin for Quantum ESPRESSO to the AiiDA platform.
The following calculations, data classes, parsers and workflows are provided:

### Calculations
- `CpCalculation`: calculation plugin for `cp.x`
- `DosCalculation`: calculation plugin for `dos.x`
- `MatdynCalculation`: calculation plugin for `matdyn.x`
- `NebCalculation`: calculation plugin for `neb.x`
- `PhCalculation`: calculation plugin for `ph.x`
- `PpCalculation`: calculation plugin for `pp.x`
- `ProjwcCalculation`: calculation plugin for `projwfc.x`
- `PwCalculation`: calculation plugin for `pw.x`
- `PwimmigrantCalculation`: to import an already completed `pw.x` into AiiDA
- `Q2rCalculation`: calculation plugin for `q2r.x`

### Data
- `ForceconstantsData`: data class for force constants produced by `q2r.x`

### Parsers
- `CpParser`: parser for the `cp.x` calculation
- `DosParser`: parser for the `dos.x` calculation
- `MatdynParser`: parser for the `matdyn.x` calculation
- `NebParser`: parser for the `neb.x` calculation
- `PhParser`: parser for the `ph.x` calculation
- `ProjwfcParser`: parser for the `projwfc.x` calculation
- `PwParser`: parser for the `pw.x` calculation
- `Q2rParser`: parser for the `q2r.x` calculation

### Workflows
- `PhBaseWorkChain`: workflow to run a `PhCalculation` to completion
- `PwBaseWorkChain`: workflow to run a `PhCalculation` to completion
- `PwRelaxWorkChain`: workflow to run a `PhCalculation` to completion
- `PwBandsWorkChain`: workflow to run a `PhCalculation` to completion
- `PwBandStructureWorkChain`: workflow to run a `PhCalculation` to completion
