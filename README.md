# AiiDA Quantum ESPRESSO
This is the official Quantum ESPRESSO plugin for AiiDA that used to be contained within the aiida_core package.
It is now maintained as a separate package and for AiiDA >= v0.9.0 the new plugin system will ensure that the
entry points are automatically registered upon installation.

# Documentation
The documentation for this package can be found on Read the Docs at 
http://aiida-quantumespresso.readthedocs.io/en/latest/

# Tests
Currently there are two test scripts included that can be used to test the basic workflows PhBaseWorkChain and PwBaseWorkChain.
You will find them in the `tests` folder. The scripts can be run directly from the commandline and they define some required
and some optional command line options. You can show these by running:

	./base.py --help

For the PwBaseWorkChain the required options are the code label of Quantum ESPRESSO's pw.x, a pseudo family name and a structure pk.
An example call would something like the following:

	./base.py -c pw-v5.4.0 -p SSSP -s 1000

If this successfully converges a PwCalculation, let's say with pk 1002, this can be used to run the PhBaseWorkChain, e.g.:

	./base.py -c ph-v5.4.0 -p 1002

# License
The aiida-quantumespresso set of plugins are released under a MIT license. See 
the LICENSE.txt file for more details.