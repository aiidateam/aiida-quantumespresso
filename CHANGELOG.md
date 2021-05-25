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
