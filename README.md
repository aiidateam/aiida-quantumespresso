# AiiDA Quantum ESPRESSO
This is the official Quantum ESPRESSO plugin for AiiDA.

The `develop` branch, which will become `aiida-quantumespresso v3.0.0`,  is only compatible with `aiida-core v1.0.0` and up.
For support for older versions of `aiida-core` use `aiida-quantumespresso v2.*`.

# Documentation
The documentation for this package can be found on Read the Docs at
http://aiida-quantumespresso.readthedocs.io/en/latest/

# Command line interface tool
The plugin comes with a builtin CLI tool: `aiida-quantumespresso`.
This tool is built using the `click` library and supports tab-completion.
To enable it, add the following to your shell loading script, e.g. the `.bashrc` or virtual environment activate script:

    eval "$(_AIIDA_QUANTUMESPRESSO_COMPLETE=source aiida-quantumespresso)"

The tool comes with various sub commands, for example to quickly launch some calculations and workchains
For example, to launch a test `PwCalculation` you can run the following command:

    aiida-quantumespresso calculation launch pw -X pw-v6.1 -p SSSP_v0.7_eff_PBE -s 134

Note that this requires the code `pw-v6.1` and pseudo potential family `SSSP_v1.1_eff_PBE` to be configured and a structure with pk `134` to be present in your database.
Each command has a fully documented command line interface, which can be printed to screen with the help flag:

    aiida-quantumespresso calculation launch ph --help

which should print something like the following:

    Usage: aiida-quantumespresso calculation launch ph [OPTIONS]

      Run a PhCalculation.

    Options:
      -X, --code CODE                 A single code identified by its ID, UUID or
                                      label.  [required]
      -C, --calculation CALCULATION   A single calculation identified by its ID or
                                      UUID.  [required]
      -k, --kpoints-mesh INTEGER...   The number of points in the kpoint mesh
                                      along each basis vector.  [default: 1, 1, 1]
      -m, --max-num-machines INTEGER  The maximum number of machines (nodes) to
                                      use for the calculations.  [default: 1]
      -w, --max-wallclock-seconds INTEGER
                                      the maximum wallclock time in seconds to set
                                      for the calculations.  [default: 1800]
      -i, --with-mpi                  Run the calculations with MPI enabled.
                                      [default: False]
      -d, --daemon                    Submit the process to the daemon instead of
                                      running it locally.  [default: False]
      -h, --help                      Show this message and exit.

# License
The aiida-quantumespresso set of plugins are released under a MIT license. See
the LICENSE.txt file for more details.

# Acknowlegements
We acknowledge support from:
* the [NCCR MARVEL](http://nccr-marvel.ch/) funded by the Swiss National Science Foundation;
* the EU Centre of Excellence "[MaX – Materials Design at the Exascale](http://www.max-centre.eu/)" (Horizon 2020 EINFRA-5, Grant No. 676598);
* the [swissuniversities P-5 project "Materials Cloud"](https://www.materialscloud.org/swissuniversities).

<img src="https://raw.githubusercontent.com/aiidateam/aiida-quantumespresso/develop/docs/source/images/MARVEL.png" width="300px" height="157px"/>

<img src="https://raw.githubusercontent.com/aiidateam/aiida-quantumespresso/develop/docs/source/images/MaX.png" width="300px" height="84px"/>

<img src="https://raw.githubusercontent.com/aiidateam/aiida-quantumespresso/develop/docs/source/images/swissuniversities.png" width="300px" height="35px"/>
