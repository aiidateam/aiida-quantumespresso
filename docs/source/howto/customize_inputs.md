# Customize inputs

## Starting from scratch

The protocols make it easier to get a starting set of inputs.
However, you can also start from an empty `builder`:

```python
pseudo_family = orm.load_group('SSSP/1.3/PBEsol/efficiency')

builder = PwBaseWorkChain.get_builder()
builder.kpoints_distance = 0.3
builder.pw.code = orm.load_code('pw@localhost')
builder.pw.structure = structure
builder.pw.pseudos = pseudo_family.get_pseudos(structure=structure)
builder.pw.parameters = {
    'SYSTEM': {
        'nbnd': 10
    },
    'CONTROL': {
        'calculation': 'scf'
    }
}
results = engine.run(builder)
```

You can also directly pass your inputs to the engine by preparing all of them in an `inputs` dictionary:

```python
pseudo_family = orm.load_group('SSSP/1.3/PBEsol/efficiency')

inputs = {
    'kpoints_distance': 0.3,
    'pw': {
        'code': orm.load_code('pw@localhost'),
        'structure': structure,
        'pseudos': pseudo_family.get_pseudos(structure=structure),
        'parameters': {
            'CONTROL': {
                'calculation': 'scf'
            }
        }
    }
}
results = engine.run(PwBaseWorkChain, **inputs)
```
