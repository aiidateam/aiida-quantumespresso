---
myst:
    substitutions:
        aiida_pseudo: '[`aiida-pseudo`](https://aiida-pseudo.readthedocs.io/)'
---

(howto-calculations-pw)=

# `pw.x`

The `pw.x` code of Quantum ESPRESSO performs many different kinds of self-consistent calculations of electronic-structure properties within Density-Functional Theory (DFT),  using a plane-wave basis set and pseudopotentials.
Examples of these properties include ground-state energy and one-electron (Kohn-Sham) orbitals, atomic forces, stresses, and structural optimization, also with variable cell.

|                     |                                                               |
|---------------------|---------------------------------------------------------------|
| Plugin class        | {class}`aiida_quantumespresso.calculations.pw.PwCalculation`  |
| Plugin entry point  | ``quantumespresso.pw``                                        |

## How to launch a `pw.x` calculation

Below is a script with a basic example of how to run a `pw.x` calculation through the {{ PwCalculation }} plugin that computes the electronic ground state of an fcc silicon crystal:

```{literalinclude} ../../tutorials/include/scripts/run_pw_basic.py
:language: python
```

Note that you may have to change the name of the code that is loaded using `load_code` and the pseudopotential family loaded with `load_group`.

## How to define input file parameters

The `pw.x` code supports many parameters that can be defined through the input file, as shown on the [official documentation](https://www.quantum-espresso.org/Doc/INPUT_PW.html).
The parameters are divided into section or "cards".
Parameters that are part of cards that start with an ampersand (`&`) should be specified through the `parameters` input of the {{ PwCalculation }} plugin.
The parameters are specified using a Python dictionary, where each card is its own sub-dictionary, for example:

```python
parameters = {
    'CONTROL': {
        'calculation': 'scf'
    },
    'SYSTEM': {
        'smearing': 'gaussian'
    },
    'ELECTRONS': {
        'electron_maxstep': 10
    }
}
```

The parameters dictionary should be wrapped in a {class}`~aiida.orm.nodes.data.dict.Dict` node and assigned to the `parameters` input of the process builder:

```python
from aiida.orm import Dict, load_code
builder = load_code('pw').get_builder()
parameters = {
    ...
}
builder.parameters = Dict(parameters)
```

:::{warning}
There are a number of input parameters that *cannot* be set, as they will be automatically set by the plugin based on other inputs, such as the `structure`.
These include:

- `CONTROL.pseudo_dir`
- `CONTROL.outdir`
- `CONTROL.prefix`
- `SYSTEM.celldm`
- `SYSTEM.nat`
- `SYSTEM.ntyp`
- `SYSTEM.a`
- `SYSTEM.b`
- `SYSTEM.c`
- `SYSTEM.cosab`
- `SYSTEM.cosac`
- `SYSTEM.cosbc`

Defining them anyway will result in an exception when launching the calculation.
:::

(howto-calculations-pw-multidimensional-parameters)=

### Multidimensional parameters

The input format of `pw.x` contains various keywords that do not simply take the format of a key value pair, but rather there will some indices in the key itself.
Take for example the [`starting_magnetization`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm287) keyword of the `SYSTEM` card.
The starting magnetization value needs to be applied to a specific species and therefore the index `i` is required to be able to make this distinction.

The {{ PwCalculation }} plugin makes this easy as it will do the conversion from kind name to species index automatically.
This allows you to specify a starting magnetization value by using a dictionary notation, where the key is the kind name to which it should be applied.
For example, if you have a structure with the kind `Co` and want it to have a given starting magnetization, one can add the following in the parameter data dictionary:

```python
parameters = {
    'SYSTEM': {
        'starting_magnetization': {
            'Co': 4.5
        }
    }
}
```

This part of the parameters dictionary will be transformed by the plugin into the following input file:

```
&SYSTEM
    starting_magnetization(1) = 4.5
/
ATOMIC_SPECIES
Co     58.93 Co.UPF
Li     6.941 Li.UPF
O      15.99 O.UPF
```

Note that since `Co` is listed as the first atomic species, the index in the `starting_magnetization(1)` keyword reflects this.
The usage of a dictionary where the keys correspond to a kind of the input structure, will work for any keyword where the index should correspond to the index of the atomic species.
Examples of keywords where this approach will work are:

- `angle1(i)`
- `angle2(i)`
- `hubbard_alpha(i)`
- `hubbard_beta(i)`
- `hubbard_j0(i)`
- `hubbard_u(i)`
- `london_c6(i)`
- `london_rvdw(i)`
- `starting_charge(i)`
- `starting_magnetization(i)`

There are also keywords that require more than index, or where the single index actually does not correspond to the index of an atomic species, such as the [`starting_ns_eigenvalue`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm466) parameters.
To allow one to define these keywords, one can use nested lists, where the first few elements constitute all the index values and the final element corresponds to the actual value.
For example the following:

```python
parameters = {
    'SYSTEM': {
        'starting_ns_eigenvalue': [
            [1, 1, 3, 3.5],
            [2, 1, 1, 2.8]
        ]
    }
}
```

will result in the following input file:

```
&SYSTEM
    starting_ns_eigenvalue(1,1,3) = 3.5
    starting_ns_eigenvalue(2,1,1) = 2.8
/
```

Note that any of the values within the lists that correspond to a kind in the input structure, will be replaced with the index of the corresponding atomic species.
For example:

```python
hubbard_j: [
    [2, 'Ni', 3.5],
    [2, 'Fe', 7.4],
]
```

would be formatted as:

```python
hubbard_j(2, 1) = 3.5
hubbard_j(2, 3) = 7.4
```

Assuming the input structure contained the kinds `Ni` and `Fe`, which would have received the atomic species indices 1 and 3 in the ultimate input file, respectively.

## How to define pseudopotentials

Each `pw.x` calculation requires a pseudopotential to be specified for each kind in the structure.
These pseudopotentials can be specified in the {{ PwCalculation }} plugin through the `pseudos` input namespace.
This input takes a dictionary, where the keys are the kind names and the values are instances of the {class}`~aiida_pseudo.data.pseudo.upf.UpfData` data plugin of the {{ aiida_pseudo }} plugin package.
For example, if the input structure is a `GaAs` crystal, the pseudopotentials could be specified as follows:

```python
from aiida.orm import load_code
from aiida_pseudo.data.pseudo.upf import UpfData

# Assume we have a UPF for Ga and As on the local file system
upf_ga = UpfData('/path/to/Ga.upf')
upf_as = UpfData('/path/to/As.upf')

builder = load_code('pw').get_builder()
builder.pseudos = {
    'Ga': upf_ga,
    'As': upf_as,
}
```

:::{tip}
We recommend using the pseudopotentials provided by the [Standard Solid-State Pseudopotentials (SSSP)](https://www.materialscloud.org/discover/sssp/table/efficiency).
The {{ aiida_pseudo }} package provides an easy and automated way to install them.
Please refer to the [section on pseudopotential setup](#installation-setup-pseudopotentials) for details.
:::

### Getting pseudopotentials from a family

If pseudopotentials have been installed as a family using the {{ aiida_pseudo }} package, they can be retrieved as follows:

```python
from ase.build import bulk
from aiida.orm import StructureData, load_code, load_group

# Load the pseudopotential family whose pseudos to use
family = load_group('SSSP/1.3/PBEsol/effiency')
structure = StructureData(ase=bulk('GaAs', 'fcc', 5.4))

builder = load_code('pw').get_builder()
builder.pseudos = family.get_pseudos_from_structure(structure=structure)
```

The {meth}`~aiida_pseudo.groups.family.pseudo.PseudoPotentialFamily.get_pseudos` method will automatically return a dictionary with the pseudos necessary for the specified `structure` that can be immediately assigned to the `pseudos` input of the builder.

### Getting recommended cutoffs from a family

Certain pseudopotential families provided by {{ aiida_pseudo }}, such as the {class}`~aiida_pseudo.groups.family.sssp.SsspFamily`, provide recommended energy cutoffs for the wave function and charge density, for each pseudopotential they provide.
Similar to the pseudopotentials themselves, these can easily be retrieved for any given `StructureData`:

```python
from ase.build import bulk
from aiida.orm import StructureData, load_code, load_group

# Load the pseudopotential family whose pseudos to use
family = load_group('SSSP/1.3/PBEsol/effiency')
structure = StructureData(ase=bulk('GaAs', 'fcc', 5.4))

builder = load_code('pw').get_builder()
cutoff_wfc, cutoff_rho = family.get_recommended_cutoffs(structure=structure, unit='Ry')
builder.parameters = {
    'SYSTEM': {
        'ecutwfc': cutoff_wfc,
        'ecutrho': cutoff_rho,
    }
}
```

Be sure to specify the `unit` as `Ry` as that is the unit that `pw.x` will expect.

## How to run an initialization-only calculation

Specify `ONLY_INITIALIZATION: True` in the `settings` input:

```python
builder = load_code('pw').get_builder()
builder.settings = Dict({'ONLY_INITIALIZATION': True})
```

If this setting is specified, the plugin will write the file `aiida.EXIT` in the working directory of the calculation.
This will cause Quantum ESPRESSO to just run the preamble of the code and then shutdown cleanly.

## How to run a gamma-only calculation

Specify `GAMMA_ONLY: True` in the `settings` input:

```python
builder = load_code('pw').get_builder()
builder.settings = Dict({'GAMMA_ONLY': True})
```

## How to fix the coordinates of atoms

Quantum ESPRESSO pw.x allows to optionally fix the coordinates of atoms during relaxation and molecular-dynamics simulations.
This functionality is enabled in {{ PwCalculation }} through the `FIXED_COORDS` setting which is a list of length equal to the number of sites in the input structure.
Each element is a list of length three containing booleans, where `True` means that the position of the site in that direction should be fixed.
For example:

```python
builder = load_code('pw').get_builder()
builder.settings = Dict({'FIXED_COORDS': [[False, False, False], [False, False, True]]})
```

will fix the position of the second site along the z-axis only.
All other coordinates are allowed to change.

## How to retrieve additional files

The {{ PwCalculation }} plugin will retrieve the most important and common output files by default.
To retrieve additional output files, specify the list of files in the `CalcJob.metadata.options.additional_retrieve_list` key in the input (**NOTE**: The usage of `ADDITIONAL_RETRIEVE_LIST` in the `settings` input is deprecated and will be removed in a future release):

```python
builder = load_code('pw').get_builder()
builder.metadata.options.additional_retrieve_list = ['custom-file.txt', 'some-other.xml']
```


## How to analyze the results

When a {{ PwCalculation }} is completed, there are quite a few possible analyses to perform.

### How to check the SCF accuracy during the self-consistent cycle

During the self-consistent field cycle, the difference in energy of the newly computed charge density and the starting one is computed and stored.
It can easily be retrieved through the {meth}`~aiida_quantumespresso.tools.calculations.pw.PwCalculationTools.get_scf_accuracy` method.
This method can be accessed directly through the `tools` of a completed {{ PwCalculation }} node:

```python
In [1]: node = load_node(IDENTIFIER)

In [2]: node.tools.get_scf_accuracy()
Out[2]:
array([1.22958881e+00, 7.84853851e-02, 5.10948147e-03, 5.42404506e-03, 2.94427169e-04, 9.25187037e-06])
```

This returns a `numpy` array with the SCF accuracy at each step during the complete cycle.
If the calculation contained dynamics, i.e. the atomic positions were updated such as in a `relax` or `vc-relax` calculation, there will be multiple SCF cycles, one for each frame of the dynamics trajectory.
To get the SCF accuracy progression of a particular frame, specify it using the `index` argument:

```python
In [1]: node = load_node(IDENTIFIER)

In [2]: node.tools.get_scf_accuracy(index=-1)
Out[2]:
array([1.38253747e+00, 5.99484465e-02, 1.20151864e-03, 4.69396365e-05, 4.08170752e-06])
```

The `index` supports negative values to start counting from the back, just as with regular `numpy` arrays.
