# AiiDA Quantum ESPRESSO
This is the official Quantum ESPRESSO plugin for AiiDA that used to be contained within the aiida_core package.
It is now maintained as a separate package and for AiiDA >= v0.9.0 the new plugin system will ensure that the
entry points are automatically registered upon installation.

# Documentation
The documentation for this package can be found on Read the Docs at 
http://aiida-quantumespresso.readthedocs.io/en/latest/

# Command line interface scripts
The plugin provides several cli scripts that make it easy to quickly launch some of the calculations and workchains
that are implemented. If you install the package with pip, the scripts will automatically be registered in your path
and you should be able to call them directly. The names can be found in the 'setup.json' under the key 'console_scripts'.

For example, to launch a test `PwCalculation` you can run the following command:

	launch_calculation_pw -c pw-v6.1 -p SSSP_v0.7_eff_PBE -s 134

Each cli script has a fully documented command line interface, which can be printed to screen with:

	launch_calculation_pw --help

which should print something like the following:

	Usage: launch_calculation_pw [OPTIONS]

	  Run a PwCalculation for a given input structure

	Options:
	  -c, --code TEXT                 the label of the AiiDA code object to use
	                                  [required]
	  -s, --structure INTEGER         the node pk of the structure  [required]
	  -p, --pseudo-family TEXT        the name of the pseudo potential family to
	                                  use  [required]
	  -k, --kpoint-mesh INTEGER...    the number of points in the kpoint mesh
	                                  along each basis vector  [default: 2, 2, 2]
	  -m, --max-num-machines INTEGER  the maximum number of machines (nodes) to
	                                  use for the calculations  [default: 1]
	  -w, --max-wallclock-seconds INTEGER
	                                  the maximum wallclock time in seconds to set
	                                  for the calculations  [default: 1800]
	  -z, --calculation-mode [scf|vc-relax]
	                                  select the calculation mode  [default: scf]
	  --help                          Show this message and exit.

# License
The aiida-quantumespresso set of plugins are released under a MIT license. See 
the LICENSE.txt file for more details.

# Acknowlegements
We acknowledge support from the [NCCR MARVEL](http://nccr-marvel.ch/) funded by the Swiss National Science Foundation and the EU Centre of Excellence "[MaX – Materials Design at the Exascale](http://www.max-centre.eu/)". (Horizon 2020 EINFRA-5, Grant No. 676598).

<img src="https://raw.githubusercontent.com/aiidateam/aiida-quantumespresso/develop/docs/source/images/MARVEL.png" width="300px" height="157px"/>

<img src="https://raw.githubusercontent.com/aiidateam/aiida-quantumespresso/develop/docs/source/images/MaX.png" width="300px" height="84px"/>
