(howto-workflows)=

# Workflows

## Default parameters

The AiiDA Quantum ESPRESSO plugin comes with a set of tested protocols. These protocols define default input parameters for different levels of precision. This concept enables the initialization of a pre-populated `builder` for most of the different `WorkChain`s, which can be obtained by calling the `get_builder_from_protocol` method.

The following example shows the working principle for the `PwBaseWorkChain`:

```
# code and structure are the mandatory input arguments

builder = PwBaseWorkChain.get_builder_from_protocol(
    code=pw_code, # load your pw-code
    structure=structure, # provide your input structure here
)

```
We refer to the documentation of the {{ get_builder_from_protocol }} method for further details about additional default input arguments.




```{toctree}
:maxdepth: 1

matdyn
ph
pw/base
pw/relax
pw/bands
q2r
pdos
```
