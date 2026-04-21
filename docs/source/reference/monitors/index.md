(reference-monitors)=

# Monitors

Monitors are callables that can be attached to a running `CalcJob` to inspect
its output files in real time and decide whether the job should be killed and
restarted (or simply stopped).  They are passed to the engine via the
`monitors` input namespace of a `BaseRestartWorkChain` subclass.

Each monitor has the signature:

```python
def monitor(node: CalcJobNode, transport: Transport) -> str | None:
    ...
```

Return `None` to let the job keep running, or a non-empty string to request
termination (the string is used as the message stored in the process report).

---

## `accuracy_stuck` — detect a plateaued accuracy metric

**Module**: `aiida_quantumespresso.tools.monitors.accuracy_stuck`

### What it does

`monitor` fetches the output file of the running calculation, scans every line
that contains the word **accuracy** (case-insensitive), and extracts the last
numerical value on that line.  If the last five parsed values are all equal
*and* consecutive, the job is considered stuck and a kill message is returned.

### Usage

```python
from aiida_quantumespresso.tools.monitors import monitor as accuracy_stuck_monitor

builder.monitors = {"accuracy_stuck": Dict({"entry_point": "quantumespresso.accuracy_stuck"})}
```

Or, when writing your own `BaseRestartWorkChain`, register the monitor via its
entry point (see `setup.cfg` / `pyproject.toml` for the registered name).

### API reference

```{eval-rst}
.. autofunction:: aiida_quantumespresso.tools.monitors.accuracy_stuck.monitor
```

```{eval-rst}
.. autofunction:: aiida_quantumespresso.tools.monitors.utils.parse_accuracy_lines
```

```{eval-rst}
.. autofunction:: aiida_quantumespresso.tools.monitors.utils.is_stuck
```
