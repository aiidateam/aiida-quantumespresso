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
numerical value on that line.  If the last `min_repeats` parsed values are all
equal *and* consecutive, the job is considered stuck and a kill message is
returned. `min_repeats` defaults to `5`.

### Usage

```python
builder.monitors = {"accuracy_stuck": Dict({"entry_point": "quantumespresso.accuracy_stuck"})}
```

To require a longer (or shorter) plateau before the job is killed, pass
`min_repeats` as a keyword argument to the monitor. As described in the
`aiida-core` how-to on
[assigning a monitor](https://aiida.readthedocs.io/projects/aiida-core/en/stable/howto/run_codes.html#how-to-assign-a-monitor),
extra keyword arguments are declared explicitly on the monitor's function
signature (never via a `**kwargs` catch-all) and are supplied through the
`kwargs` key of the settings `Dict`:

```python
builder.monitors = {
    "accuracy_stuck": Dict({
        "entry_point": "quantumespresso.accuracy_stuck",
        "kwargs": {"min_repeats": 10},
    })
}
```

Or, when writing your own `BaseRestartWorkChain`, register the monitor via its
entry point (see `setup.cfg` / `pyproject.toml` for the registered name).
