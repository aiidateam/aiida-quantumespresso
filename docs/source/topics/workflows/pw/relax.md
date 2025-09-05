(topics-workflows-pw-relax)=

# `PwRelaxWorkChain`

**Purpose:** Obtain the geometric ground state of a structure.

## Logic

The `PwRelaxWorkChain` consist of two main components:

1. An initial relaxation with looser precision settings.
2. A full relaxation loop.

The initial relaxation makes it less costly to obtain a reasonable first ground state geometry, improving the efficiency and the less strict parameters makes it easier to converge for geometries far removed from the ground state, improving robustness.

If the initial geometry optimization completes successfully, the optimization workflow is run again but this time with the convergence parameters as determined by the input protocol.

The full relaxation step is repeated until the following two conditions are met:

* **No Pulay stress:** The stresses in the final SCF performed by `pw.x` are below the selected convergence threshold.
  During the geometry optimization, the basis sets used (i.e. the list of reciprocal-space G vectors), that depend not only on the energy cutoffs but also on the crystal geometry, are not updated at each step.
  The final SCF step that is performed at the end of the optimization recomputes the basis sets on the optimized unit cell.
  Large stresses in the final SCF indicate that the optimized geometry differs significantly from the initial geometry and the static basis sets used towards the end of the optimization were no longer converged enough.
* **Sufficient k-mesh density:** The k-point mesh corresponds to a lower density than the one dictated by the protocol.
  Since the unit cell can change during the optimization, it is possible that, towards the end of the optimization cycle, the initial k-point mesh no longer satisfies the minimum required k-point density specified by the input protocol.
  In this case, a new k-point mesh is generated and another geometry optimization is performed.

## Input schema

```{eval-rst}
.. aiida-workchain:: PwRelaxWorkChain
    :module: aiida_quantumespresso.workflows.pw.relax
```
