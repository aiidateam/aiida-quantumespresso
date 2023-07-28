---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: QE-dev
    language: python
    name: qe-dev
---

(tutorials-magnetic-configurations)=

# Magnetic configurations

There are several ways of defining a magnetic configuration in Quantum ESPRESSO.
This tutorial focuses on using {{ starting_magnetization }} to specify the initial magnetization for each kind in our system.

(tutorials-magnetic-configuration-parameters)=

## Using the `parameters`

One way to specify the initial magnetic configuration of your system is by simply using the `parameters` input of the {{ PwCalculation }}.
Start by loading your default profile:

```python
from aiida import orm, load_profile

load_profile()
```

Next, set up a basic structure for HCP cobalt, and load the code you have set up for `pw.x`:

:::{margin}
❗️Make sure to replace the label of the `pw.x` code in the {func}`~aiida.orm.load_code` function by that of _your_ installed code.
:::

```python
from ase.build import bulk

structure = orm.StructureData(ase=bulk('Co', 'hcp', 2.5, 4.06))
pw_code = orm.load_code('qe-7.2-pw@localhost')
```

Obtain a pre-populated builder using the {{ get_builder_from_protocol }} method:

:::{margin}
</br></br></br>
🚀 Use the `fast` protocol for testing and demonstrations!
:::

```python
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

builder = PwBaseWorkChain.get_builder_from_protocol(
    structure=structure,
    code=pw_code,
    protocol='fast'
)
```

In order to set the magnetic moments, you need to change the `parameters` input of the {{ PwCalculation }}:

```python
parameters = builder.pw.parameters.get_dict()
parameters['SYSTEM']['nspin'] = 2
parameters['SYSTEM']['starting_magnetization'] = {'Co': 1}
builder.pw.parameters = orm.Dict(parameters)
```

The commands above do the following:

1.  Get the `parameters` input of the `pw` namespace as a regular Python `dict`.
2.  Set {{ nspin }} to `2` so Quantum ESPRESSO does a spin-polarised calculation.
3.  Set the {{ starting_magnetization }} for the `Co` atomic species or kind to 1.
4.  Replace the `parameters` of the {{ PwCalculation }} by a new `Dict` node that has the {{ starting_magnetization }} set.

:::{tip}
Note that {{ starting_magnetization }} is a "multidimensional parameter", which in the `aiida-quantumespresso` plugin you shouldn't specify using indices between brackets as you would in a regular Quantum ESPRESSO input file.
See the [corresponding how-to section](howto-calculations-pw-multidimensional-parameters) for more details.
:::

Now we ready to submit the builder to the AiiDA daemon:

:::{margin}
⏱ Since this cell actually runs a `pw.x` calculation, it might take a bit longer to complete.
:::

```python
from aiida.engine import run_get_node

results, pwbase_node = run_get_node(builder)
```

Once the calculation is completed, you can for example get the total magnetization in the unit cell from the output parameters:

```python
results['output_parameters'].get_dict()['total_magnetization']
```

:::{important}
The {{ starting_magnetization }} input is defined as "spin polarization" of the atom species, with inputs ranging from -1 (all spins down) to 1 (all spins up).
However, the "magnetization" outputs are always defined in units of Bohr magneton.
We will commonly refer to the values in Bohr magneton as the "magnetic moments" to differentiate them from "magnetization".
:::

(tutorials-magnetic-configuration-spintype)=

## Using the {{ SpinType }}

An alternative way of performing spin-polarised calculations is by specifying a {{ SpinType }} and passing it to the {{ get_builder_from_protocol }} method:

```python
from aiida_quantumespresso.common.types import SpinType

builder = PwBaseWorkChain.get_builder_from_protocol(
    structure=structure,
    code=pw_code,
    protocol='fast',
    spin_type=SpinType.COLLINEAR
)
```

In this case, an initial {{ starting_magnetization }} will be set for each element based on the following logic:

-   If the element has a partially occupied *d* or *f* shell, calculate what the {{ starting_magnetization }} should be in case we aim for a maximally possible magnetic moment, i.e. 5 or 7 Bohr magneton for *d* or *f* shells.
-   For other elements, simply set {{ starting_magnetization }} to 0.1.


:::{note}
The approach outlined above will always set *positive* values for {{ starting_magnetization }}, and hence can only initialise the system in a ferro or ferrimagnetic state.
:::

You can have a look to see what {{ starting_magnetization }} was set by inspecting the contents of the `parameters` set in the builder:

```python
builder.pw.parameters.get_dict()['SYSTEM']['starting_magnetization']
```

Once again, try running the calculation:

```python
results, pwbase_node = run_get_node(builder)
```
The calculation should have the same final total magnetization as the one in [section using the parameters](tutorials-magnetic-configuration-parameters).

```python
results['output_parameters'].get_dict()['total_magnetization']
```

## Anti-ferromagnetic configurations


What if you want to initialize the system in an anti-ferromagnetic state?
In this case, Quantum ESPRESSO requires two different kinds to be specified in the inputs, but looking at the `sites` in the input structure:

```python
structure.sites
```

You can see that although our structure has two _sites_, they are both of the same _kind_.
This means you'll have to create a new `` node with two _different_ kinds.
For this purpose, you can use the {{ create_magnetic_configuration }} calculation function:

```python
from aiida_quantumespresso.calculations.functions.create_magnetic_configuration import create_magnetic_configuration

results = create_magnetic_configuration(structure, [-2, 2])
```

Note that as {{ create_magnetic_configuration }} is a calculation function, and so its execution is tracked in the provenance:

```console
❯ verdi process list -ap1
  PK  Created    Process label                 Process State    Process status
----  ---------  ----------------------------  ---------------  ----------------
 [...]
 695  12s ago    create_magnetic_configuration  ⏹ Finished [0]

Total results: 8

Report: last time an entry changed state: 3m ago (at 02:06:20 on 2023-05-28)
```

You can see the in- and outputs of the {{ create_magnetic_configuration }} process using `verdi process show`:

```bash
❯ verdi process show 695
Property     Value
-----------  ------------------------------------
type         create_magnetic_configuration
state        Finished [0]
pk           695
uuid         bb8d4026-19c2-4efd-a99d-b55336df0eb3
label        create_magnetic_configuration
description
ctime        2023-05-28 04:50:25.887472+02:00
mtime        2023-05-28 04:50:26.098163+02:00

Inputs                      PK  Type
------------------------  ----  -------------
atol                       693  Float
magnetic_moment_per_site   692  List
structure                  659  StructureData
ztol                       694  Float

Outputs             PK  Type
----------------  ----  -------------
magnetic_moments   697  Dict
structure          696  StructureData
```

:::{important}
Every time you execute the {{ create_magnetic_configuration }} calculation function, it will store another process with corresponding in- and output nodes to the database!
If you want to avoid this, you can set `store_provenance` to `False` in the `metadata`:

```python
results = create_magnetic_configuration(structure, [-2, 2], metadata={"store_provenance": False})
```
:::

Let's have a look at the `structure` output:

```python
new_structure = results['structure']
new_structure
```

As you can see, this is the new {{ StructureData }}, which this time has two different kind names for each of the `Co` sites:

```python user_expressions=[]
results['structure'].sites
```

The `results` also contain the `magnetic_moments` output:

```python
results['magnetic_moments'].get_dict()
```

This can be passed to the `initial_magnetic_moments` input of the {{ get_builder_from_protocol }} method, along with the `new_structure`:

```python
builder = PwBaseWorkChain.get_builder_from_protocol(
    structure=new_structure,
    code=pw_code,
    protocol='fast',
    spin_type=SpinType.COLLINEAR,
    initial_magnetic_moments=results['magnetic_moments']
)
```

Let's have another look at the {{ starting_magnetization }} input:

```python
builder.pw.parameters['SYSTEM']['starting_magnetization']
```

As always, we can run the calculation:

```python
results, pwbase_node = run_get_node(builder)
```

In this case, the _total_ magnetization in the unit cell is zero:

```python
results['output_parameters'].get_dict()['total_magnetization']
```

But the _absolute_ magnetization is not:

```python
results['output_parameters'].get_dict()['absolute_magnetization']
```

## Passing the magnetic configuration

Often you'll want to pass the magnetic configuration obtained from one calculation to the next.
The final magnetic configuration of the `pw.x` run can be found in the `output_trajectory` output:

```python
trajectory = pwbase_node.outputs.output_trajectory
```

This {class}`aiida.orm.TrajectoryData` node contains a lot of information for each _ionic_ step performed in the `pw.x` calculation, all stored as `numpy` arrays in the repository.
To see all the arrays stored you can use the `get_arraynames()` method:

```python
trajectory.get_arraynames()
```

We're interested in the `atomic_magnetic_moments`.
Retrieve the _final_ magnetic moments using the `get_array()` method:

```python
magnetic_moments = trajectory.get_array('atomic_magnetic_moments')[-1]
```

:::{note}
Here you retrieved the _final_ array of magnetic moments by specifying the `-1` index.
Of course, for a static (`scf`) calculation there is only one ionic step in the trajectory.
If you had run a `vc-relax` calculation, the arrays in `output_trajectory` would contain all the information above for _each_ ionic step.
:::

Have a look at the contents of `magnetic_moments`.
In this case, the magnetic moments are stored as a one-dimensional array of length 2.
Note that Quantum ESPRESSO prints the magnetic configuration as a list of magnetic moments:

```fortran
     Magnetic moment per site  (integrated on atomic sphere of radius R)
     atom   1 (R=0.411)  charge= 15.6004  magn=  1.2090
     atom   2 (R=0.411)  charge= 15.6004  magn= -1.2087
```

And that, same as for the total magnetization, the magnetic moments are expressed in Bohr magneton.

Due to this discrepancy between the in- and outputs in Quantum ESPRESSO, it becomes a bit challenging to properly pass magnetic configurations between calculations.
You'd have to perform roughly the following steps:

1. Retrieve the final structure and atomic moments.
   Which structure to use will depend on the calculation you have run.
2. Check if the structure needs to have different kinds defined in order to specify the magnetic moments.
3. Properly map the kinds to the the correct magnetic moments, grouping together sites whose magnetic moments are close.
3. Check if all magnetic moments are very small.
   In this case continue with a non-magnetic calculation.

Fortunately, the {{ PwCalculation }} comes with a handy tool for obtaining the magnetic configuration.
First, obtain the final {{ PwCalculation }} from the `called` processes of the {{ PwBaseWorkChain }} node:

```python
pw_node = pwbase_node.called[-1]
```

Then, you can simply use the {{ get_magnetic_configuration }} method in the `tools` name space:

```python
configuration = pw_node.tools.get_magnetic_configuration()
```

:::{note}
In case new kinds are needed for the configuration, the {{ create_magnetic_configuration }} calculation function is run inside {{ get_magnetic_configuration }} to generate the new structure.
This step is added to the provenance and hence stored in the database.
:::

The configuration contains both the new `structure` and `magnetic_moments`, which can then be easily passed to the {{ get_builder_from_protocol }} method:

```python
builder = PwBaseWorkChain.get_builder_from_protocol(
    structure=configuration.structure,
    code=pw_code,
    spin_type=SpinType.COLLINEAR,
    initial_magnetic_moments=configuration.magnetic_moments
)
```

Have a look at the {{ starting_magnetization }} input that was generated:

```python
builder.pw.parameters['SYSTEM']['starting_magnetization']
```

Is it different from the one at the start of the initial run?
