(howto-workflows-pw-base)=

# `PwBaseWorkChain`

To obtain a builder with the default protocol inputs:

```
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

builder = PwBaseWorkChain.get_builder_from_protocol(
    code=orm.load_code('pw@localhost'),
    structure=orm.load_node(72894),
)
```

The calculation can then be run using for example the `run` function of the engine:

```
from aiida import engine

engine.run(builder)
```
