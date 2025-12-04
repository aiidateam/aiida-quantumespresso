# Workflow logic

## `PwBaseWorkChain`
**Purpose:** Run a `pw.x` calculation with error handling.

The `PwBaseWorkChain` provides automated error handling for common failures encountered during Quantum ESPRESSO `pw.x` calculations.
When an error is detected, the work chain attempts to fix the issue and restart the calculation automatically.
The following sections describe each error handler and the recovery strategy employed.

### Band occupation sanity check

- **Handler:** `sanity_check_insufficient_bands`
- **Exit codes handled:** Successful calculations (exit code 0)

This handler performs a post-calculation check on successfully converged calculations to check if for any spin channel or k-point the highest band has an occupation > 0.005, which indicates that insufficient bands have been used in the calculation.
In this case, increase the number of bands (`nbnd`) by 5% with a minimum of 4 bands, and restart from the charge density.

### Diagonalization errors

- **Handler:** `handle_diagonalization_errors`
- **Exit codes handled:**
  - `ERROR_COMPUTING_CHOLESKY` (462)
  - `ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED` (463)
  - `ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE` (464)
  - `ERROR_ZHEGVD_FAILED` (465)
  - `ERROR_QR_FAILED` (466)
  - `ERROR_EIGENVECTOR_CONVERGENCE` (467)
  - `ERROR_BROYDEN_FACTORIZATION` (468)

When the diagonalization algorithm fails, this handler systematically tries alternative algorithms in order of robustness.
The default algorithm is `david`, and alternatives are tried in the order: `ppcg`, `paro`, and finally `cg` (conjugate gradient).
If all algorithms have been tried, the calculation aborts.

### Out of walltime

- **Handler:** `handle_out_of_walltime`
- **Exit codes handled:**
  - `ERROR_OUT_OF_WALLTIME` (400)

Handles calculations that exceeded the allocated walltime but shut down cleanly.
Uses the output structure if available (for relaxation calculations) and perform a full restart with `CONTROL.restart_mode = 'restart'`.

### Ionic cycle interruption

- **Handler:** `handle_ionic_interrupted_partial_trajectory`
- **Exit codes handled:**
  - `ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY` (503)

Handles calculations that were interrupted during ionic optimization due to transient problems, where the charge density and wave functions are likely corrupt.
Uses the last output structure and restarts from scratch.

### Ionic convergence reached except in final SCF

- **Handler:** `handle_vcrelax_converged_except_final_scf`
- **Exit codes handled:**
  - `ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF` (501)

For `vc-relax` calculations where ionic convergence thresholds are met but the final SCF exceeds thresholds, marks the workchain as finished and returns the output structure with an exit code to indicate partial success.
This information is then used by the `PwRelaxWorkChain`, whose job is to remove all Pulay stresses from the system.

### BFGS history failures

- **Handler:** `handle_relax_recoverable_ionic_convergence_bfgs_history_error`
- **Exit codes handled:**
  - `ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE` (520)
  - `ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE` (521)

Handles failures in the BFGS ionic minimization where the history is reset twice consecutively.
For `relax` calculations, switches to `damp` dynamics; for `vc-relax`, either reduces `trust_radius_min` by 90% (if above 1.0e-4) or switches to `damp`/`damp-w` dynamics.
All cases restart from the previous charge density using the output structure.

### Ionic convergence failures

- **Handler:** `handle_relax_recoverable_ionic_convergence_error`
- **Exit codes handled:**
  - `ERROR_IONIC_CONVERGENCE_NOT_REACHED` (500)
  - `ERROR_IONIC_CYCLE_EXCEEDED_NSTEP` (502)
  - `ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE` (520, fallback)
  - `ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE` (521, fallback)

Handles cases where ionic convergence thresholds were not met but the calculation shut down cleanly with a usable output structure.
Restarts from the previous charge density using the last output structure.

### Significant volume contraction in vc-relax

- **Handler:** `handle_vcrelax_recoverable_fft_significant_volume_contraction_error`
- **Exit codes handled:**
  - `ERROR_RADIAL_FFT_SIGNIFICANT_VOLUME_CONTRACTION` (542)

Handles significant volume scaling during variable-cell relaxation, which requires recalculation of pseudopotential tables.
Uses the output structure, doubles the `CELL.cell_factor` parameter, and restarts from scratch.

### Electronic convergence failures in relaxations

- **Handler:** `handle_relax_recoverable_electronic_convergence_error`
- **Exit codes handled:**
  - `ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED` (510)
  - `ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED` (511)

Handles electronic convergence failures during ionic relaxation where the output structure is still usable.
Reduces `mixing_beta` by 20%, uses the output structure, and restarts from scratch.

### Electronic convergence not reached

- **Handler:** `handle_electronic_convergence_not_reached`
- **Exit codes handled:**
  - `ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED` (410)

Handles electronic convergence failures in single-point or band structure calculations.
Reduces `mixing_beta` by 20% and performs a full restart from the previous calculation.

### Electronic convergence warnings

- **Handler:** `handle_electronic_convergence_warning`
- **Exit codes handled:**
  - `WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED` (710)

Handles cases where electronic convergence was not reached but the input parameters indicate this is acceptable (e.g., `scf_must_converge = False` or `electron_maxstep = 0`).
Marks the workchain as finished and returns the outputs with a warning exit code.


## `PwRelaxWorkChain`
**Purpose:** Obtain the geometric ground state of a structure.

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
