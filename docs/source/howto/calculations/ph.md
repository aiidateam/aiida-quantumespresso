(howto-calculations-ph)=

# `ph.x`

The `ph.x` code of Quantum ESPRESSO is used to compute phonons using density-functional perturbation theory.

|                     |                                                               |
|---------------------|---------------------------------------------------------------|
| Plugin class        | {class}`aiida_quantumespresso.calculations.ph.PhCalculation`  |
| Plugin entry point  | ``quantumespresso.ph``                                        |

## How to launch a `ph.x` calculation

:::{note}
In order to run a `ph.x` calculation, you first need to have completed a `pw.x` calculation.
See the [tutorial](#tutorials-pw-through-api) or [how-to guide](#howto-calculations-pw) for more information.
:::

Once you have successfully run a `PwCalculation` you can run a `ph.x` calculation through the `PhCalculation` plugin as follows:

```{literalinclude} ../../tutorials/include/scripts/run_ph_basic.py
:language: python
```

Note that you will have to replace `IDENTIFIER_PW_CALCULATION` with the identifier (pk or UUID) of the completed `PwCalculation`.

## How to define input file parameters

The `ph.x` code supports many parameters that can be defined through the input file, as shown on the [official documentation](https://www.quantum-espresso.org/Doc/INPUT_PH.html).
Parameters that are part of the `INPUTPH` card should be specified through the `parameters` input of the `PwCalculation` plugin.
The parameters are specified using a Python dictionary, for example:

```python
parameters = {
    'INPUTPH': {
        'tr2_ph' : 1.0e-8,
        'epsil' : True,
        'ldisp' : True,
    }
}
```

The parameters dictionary should be wrapped in a {py:class}`~aiida.orm.nodes.data.dict.Dict` node and assigned to the `parameters` input of the process builder:

```python
from aiida.orm import Dict, load_code
builder = load_code('ph').get_builder()
parameters = {
    ...
}
builder.parameters = Dict(parameters)
```

The q-points of the input file are specified with a `KpointsData` node through the `qpoints` input of the `PhCalculation` plugin.

:::{warning}
There are a number of input parameters that *cannot* be set, as they will be automatically set by the plugin based on other inputs, such as the `structure`.
These include:

- `INPUTPH.outdir`
- `INPUTPH.verbosity`
- `INPUTPH.prefix`
- `INPUTPH.fildyn`
- `INPUTPH.ldisp`
- `INPUTPH.nq1`
- `INPUTPH.nq2`
- `INPUTPH.nq3`
- `INPUTPH.qplot`
:::
