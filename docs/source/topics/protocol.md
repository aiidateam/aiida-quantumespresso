# Protocol

In the context of this package, a _protocol_ is a set of input parameters used to run a workflow.
You can obtain a fully populated builder for each `WorkChain` by providing one of the three supported protocols to the `get_builder_from_protocol()` method:

- `fast`: low precision calculations at minimal computational cost for testing purposes.
- `balanced`: normal precision calculations at moderate computational cost.
- `stringent`: high precision calculations at higher computational cost.
  Recommended for metals with lanthanides/actinides

On this page you can find an overview of the various input parameters that are controlled by protocols and their values.

## K-points sampling

The Brillouin zone is sampled at $k$-points that are defined by a Monkhorst-Pack mesh including the $\Gamma$-point.
The mesh density is defined in terms of the `kpoints_distance` input, which defined the maximum distance between $k$-points in each reciprocal-space direction (i.e., the protocol chooses the smallest $k$-point mesh with at least the density required by the specified `kpoints_distance`).
These values correspond to the extensively tested `balanced` protocol described in detail by Nascimento _et al._[^gabby].

[^gabby]: Nascimento, G., dos Santos, F. J., Bercx, M., Grassano, D., Pizzi, G., & Marzari, N., [Accurate and efficient protocols for high-throughput first-principles materials simulations](https://arxiv.org/abs/2504.03962), _preprint_ (2025).

| Protocol name | Smearing [Ry]   | `kpoints_distance` [1/A] |
|---------------|:---------------:|:-----------------------:|
| `fast`        | 0.0275          | 0.30                    |
| `balanced`    | 0.0200          | 0.15                    |
| `stringent`   | 0.0125          | 0.10                    |

## Pseudopotentials

Pseudopotentials are taken from the [Standard Solid-state Pseudopotential (SSSP)](https://www.materialscloud.org/discover/sssp/table/efficiency) library, which collects pseudopotentials from a number of libraries.
The SSSP provides a set of rigorously tested values for the recommended wave function and charge density energy cutoffs for each pseudopotential.
For every structure, [`ecutwfc`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id51) and [`ecutrho`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id52) are set to the maximum of the SSSP-recommended values over the elements in its composition.

| Protocol name | Pseudos & Cutoffs            |
|---------------|------------------------------|
| `fast`        | [`SSSP/1.3/PBEsol/efficiency`](https://www.materialscloud.org/discover/sssp/table/efficiency)  |
| `balanced`    | [`SSSP/1.3/PBEsol/efficiency`](https://www.materialscloud.org/discover/sssp/table/efficiency)  |
| `stringent`   | [`SSSP/1.3/PBEsol/precision`](https://www.materialscloud.org/discover/sssp/table/precision)    |

```{note}
As by default in Quantum ESPRESSO the exchange-correlation functional is taken from the pseudopotential files, this also means all protocols use PBEsol.
```

## Thresholds

The thresholds for electronic and ionic convergence are the following:

- SCF energy: the energy threshold for self-consistency in the SCF cycle, see [`conv_thr`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id149).
  The protocol is defined in terms of Ry/atom, i.e. the value of `conv_thr` is calculated by multiplying the protocol value with the number of atoms.
- Ionic energy: the energy threshold for convergence in the ionic optimization, see [`etot_conv_thr`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id18).
  The protocol is defined in terms of Ry/atom, i.e. the value of `etot_conv_thr` is calculated by multiplying the protocol value with the number of atoms.
- Forces: the threshold for convergence of the force components, see [`forc_conv_thr`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id19).

| Protocol name | SCF energy [Ry/atom] | Ionic energy [Ry/atom]  | Forces [Ry/bohr]  |
|---------------|:--------------------:|:-----------------------:|:-----------------:|
| `fast`        | 4e-10                | 1e-4                    | 1e-3              |
| `balanced`    | 2e-10                | 1e-5                    | 1e-4              |
| `stringent`   | 1e-10                | 5e-6                    | 5e-5              |

## Magnetism 

```{note}
See the [tutorial on magnetic calculations](tutorials-magnetic-configurations) for an introduction on how to work with magnetism in `aiida-quantumespresso`.
```

The `get_builder_from_protocol()` method supports defining the type of magnetic calculation via the `spin_type` input:

- `SpinType.NONE`: Non-spin-polarised calculation ([`nspin`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id73) = 1).
- `SpinType.COLLINEAR`: Spin-polarised calculation ([`nspin`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id73) = 2).
  Each kind is initialized in a high-spin ferromagnetic configuration, where elements with partially occupied $d$ or $f$ orbitals are assigned a magnetic moment of 5 $\mu_B$ or 7 $\mu_B$, respectively.
  For all other elements, the electrons are initialized to have a 10% surplus in the spin-up channel ([`starting_magnetization`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id50) = 0.1).
- `SpinType.NON_COLLINEAR`: Non-collinear spin-polarised calculation ([`noncolin`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id79) = `True`, [`nspin`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id73) = 4).
  The size of the magnetic vector ([`starting_magnetization`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id50)) is set to the same value as the collinear case, and both angles ([`angle1`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id104) and [`angle2`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id105)) are set to zero. 
  We also use the fully relativistic PBEsol pseudopotentials from [the Pseudo Dojo](https://www.pseudo-dojo.org/) (`PseudoDojo/0.4/PBEsol/FR/standard/upf` family).
- `SpinType.SPIN_ORBIT`: Non-collinear spin-polarised calculation including spin-orbit [`lspinorb`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id111).
  The ([`starting_magnetization`](https://www.quantum-espresso.org/Doc/INPUT_PW.html#id50)) and angles are initialised using the same approach as for `SpinType.NON_COLLINEAR`.
  