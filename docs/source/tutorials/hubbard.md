(tutorials-hubbard)=

# DFT with Hubbard corrections

Since Quantum ESPRESSO v7.1, a new syntax for Hubbard corrections has been introduced.

::: {note}
Hubbard parameters should be computed from first-principles, e.g. using linear-response theory.
You can use the {{ aiida_hubbard }} package to calculate $U$ and $V$ parameters self-consistently from first-principles.

Please, refer also to our last publication and cite it if you make use of the {{ HubbardStructureData }} in any part of your work:

> Lorenzo Bastonero _et al._, [_First-principles Hubbard parameters with automated and reproducible workflows_](https://arxiv.org/abs/2503.01590), arXiv:2503.01590v1 (2025)
:::


## Defining the Hubbard corrections through the `HubbardStructureData`

The Hubbard correction is a corrective term that is added to the Hamiltonian of a system
which suffers from great __self-interaction errors__. This is usually the case for transition
metals on their _d_ manifolds. An extra correction to account for the hybridization can be accounted
for with the ligands, typically belonging to the _p_ element group. Such interaction needs to be
localized in space. This is the reason why we need to define the __projectors__. Quantum ESPRESSO
allows you to define different type of projections $| \phi^I_m \rangle$ ($m$ orbital quantum number, $I$ atom in cell). Finally, one can work with different formulations (e.g. Dudarev or Liechtenstein).

Therefore, for a Hubbard corrected calculation we need to define:

1. _Which atoms_ $I$ and $J$ and _which manifolds_ $l$ and $l'$, defining $V^{IJ}$ ($U^I \equiv V^{II}$).
2. The kind of _projectors_ (default to _ortho-atomic_).
3. The formulation (default to _Dudarev_).

Since manifolds and interacting atoms belong to the structure, then you need to definet them together as an {{ HubbardStructureData }}.

In the following, we take LiCoO{sub}`2` as example, and we suppose we want to target the _3d_ orbitals of cobalt and the intersite interaction between _2p_ of oxygen and _3d_ of cobalt.


```python
# Load the AiiDA profile
from aiida import orm, load_profile
load_profile()
```




    Profile<uuid='bb03446a75234a54b23c64546f6c8712' name='default'>



Let's define the {{ HubbardStructureData }}:

:::{note}
:class: dropdown

If you already have a {{ StructureData }}, you can load the structure information in {{ HubbardStructureData }} as follows:

```python
my_structure = load_node(IDENTIFIER)
hubbard_structure = HubbardStructureData.from_structure(my_structure)
```
:::


```python
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

a, b, c, d = 1.40803, 0.81293, 4.68453, 1.62585
cell = [[a, -b, c], [0.0, d, c], [-a, -b, c]]
sites = [
    # symbol, "kind name", position
    ['Co', 'Co', (0, 0, 0)],
    ['O',   'O', (0, 0, 3.6608)],
    ['O',   'O', (0, 0, 10.392)],
    ['Li', 'Li', (0, 0, 7.0268)],
]
hubbard_structure = HubbardStructureData(cell=cell, sites=sites)
hubbard_structure.initialize_onsites_hubbard(
    atom_name="Co",
    atom_manifold="3d",
    value=5.0,
)
hubbard_structure.initialize_intersites_hubbard(
    atom_name="Co",
    atom_manifold="3d",
    neighbour_name="O",
    neighbour_manifold="2p",
    value=1.0e-6,
)
```

Let's visualize what will be print in the Hubbard card of Quantum ESPRESSO.


```python
from aiida_quantumespresso.utils.hubbard import HubbardUtils

print(HubbardUtils(hubbard_structure).get_hubbard_card())
```

    HUBBARD	ortho-atomic
     V	Co-3d	Co-3d	1	1	5.0
     V	Co-3d	O-2p	1	46	1e-06



As you can see, the desired interactions have been initialized correctly.

:::{important}
If you need to run the `hp.x` binary, you need to make sure to have your 'Hubbard atoms' first in the list of atoms. This is due to the way the `hp.x` routine works internally, requiring those to be first.
You can simply do this with the following snippet (IF THE NODE IS NOT YET STORED!):

```python
from aiida_quantumespresso.utils.hubbard import HubbardUtils
HubbardUtils(hubbard_structure).reorder_atoms
```
:::

::: {note}
The {class}`~aiida_hubbard.workflows.hubbard.SelfConsistentHubbardWorkChain` of the {{ aiida_hubbard }} plugin takes care of calculating $U$ and $V$ parameters from first-priciples are accounting for all this internal Quantum ESPRESSO details (including parallelization over atoms, **q**-points, renaming of the atoms, etc).
:::

### Modify projectors and formulation

Here is how you can modify the Hubbard projectors and formulation.


```python
hubbard = hubbard_structure.hubbard
hubbard.projectors = 'atomic' # default is 'ortho-atomic'
hubbard.formulation = 'dudarev' # already the default
hubbard_structure.hubbard = hubbard
print(HubbardUtils(hubbard_structure).get_hubbard_card())
```

    HUBBARD	atomic
     V	Co-3d	Co-3d	1	1	5.0
     V	Co-3d	O-2p	1	46	1e-06




```python
hubbard = hubbard_structure.hubbard
hubbard.projectors = 'ortho-atomic'
hubbard_structure.hubbard = hubbard
print(HubbardUtils(hubbard_structure).get_hubbard_card())
```

    HUBBARD	ortho-atomic
     V	Co-3d	Co-3d	1	1	5.0
     V	Co-3d	O-2p	1	46	1e-06



## Other ways to define the {{ HubbardStructureData }}

We show here two other ways to define the Hubbard parameters.

### Define the Hubbard $U$ and $V$ using nearest neighbour analysis


```python
from aiida_quantumespresso.utils.hubbard import initialize_hubbard_parameters

hubbard_structure = initialize_hubbard_parameters(
    structure=hubbard_structure, # we use the previous hubbard_structure for convenience, but can be a StructureData
    pairs={'Co': ['3d', 5.0, 1.0, {'O':'2p'}]},
    fold=True, # whether to fold the atoms in the unit cell (convenient)
)

print(HubbardUtils(hubbard_structure).get_hubbard_card())
```

    HUBBARD	ortho-atomic
     V	Co-3d	Co-3d	1	1	5.0
     V	Co-3d	O-2p	1	11	1.0
     V	Co-3d	O-2p	1	43	1.0
     V	Co-3d	O-2p	1	19	1.0
     V	Co-3d	O-2p	1	46	1.0
     V	Co-3d	O-2p	1	22	1.0
     V	Co-3d	O-2p	1	54	1.0



### Parse Hubbard parameters from file

The `hp.x` code would return a file called `HUBBARD.dat` that you can parse inside the structure.

::: {note}
If you have a `pw.x` input file, you can also parse the file to get the structure using the following snippet

```python
from aiida_quantumespresso.tools.pwinputparser import PwInputFile

with open('pw.in', 'r') as handle:
    parser = PwInputFile(handle.read())
structure = parser.get_structuredata()
```
:::


```python
import tempfile

text = """HUBBARD	ortho-atomic
 V	Co-3d	Co-3d	1	1	5.0
 V	Co-3d	O-2p	1	11	1.5
 V	Co-3d	O-2p	1	43	1.5
 V	Co-3d	O-2p	1	19	1.5
 V	Co-3d	O-2p	1	46	1.5
 V	Co-3d	O-2p	1	22	1.5
 V	Co-3d	O-2p	1	54	1.5
"""

hubbard_structure.clear_hubbard_parameters() # remove all Hubbard parameters
hubbard_utils = HubbardUtils(hubbard_structure)

with tempfile.NamedTemporaryFile(mode="w+") as handle:
    handle.write(text)
    handle.seek(0)
    hubbard_utils.parse_hubbard_dat(handle.name)

print(hubbard_utils.get_hubbard_card())
```

    HUBBARD	ortho-atomic
     V	Co-3d	Co-3d	1	1	5.0
     V	Co-3d	O-2p	1	11	1.5
     V	Co-3d	O-2p	1	43	1.5
     V	Co-3d	O-2p	1	19	1.5
     V	Co-3d	O-2p	1	46	1.5
     V	Co-3d	O-2p	1	22	1.5
     V	Co-3d	O-2p	1	54	1.5



Here we used a temporary file, but you can simply change the code to

```python
hubbard_utils.parse_hubbard_dat('/path/to/HUBBARD.dat')
```

## Run a $DFT+U+V$ calculation

Now that we defined a {{ HubbardStructureData }}, we can use it as a {{ StructureData }} _anywhere_ in the `aiida-quantumespresso` world (e.g. also in other packages such as [`aiida-vibroscopy`](https://aiida-vibroscopy.readthedocs.io/en/latest/5_iraman_functionals.html) to compute phonons, IR and Raman spectra).

::: {attention}
Make sure to:
* Adapt the `pw_code`
* Have at least Quantum ESPRESSO v7.1
:::


```python
from aiida.engine import run_get_node
from aiida_quantumespresso.common.types import SpinType
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

# !!! Make sure to adapt the `pw_code` !!!
# !!! You need at least QE v7.1        !!!
pw_code = 'pw@localhost'
builder = PwBaseWorkChain.get_builder_from_protocol(
    structure=hubbard_structure,
    code=pw_code,
    protocol='fast',
    spin_type=SpinType.COLLINEAR, # usually transition metals need spin polarized calculations
)

results, calc = run_get_node(builder)
```

    05/09/2025 09:27:13 AM <99094> aiida.orm.nodes.process.workflow.workchain.WorkChainNode: [REPORT] [106254|PwBaseWorkChain|run_process]: launching PwCalculation<106259> iteration #1
    05/09/2025 09:42:47 AM <99094> aiida.orm.nodes.process.workflow.workchain.WorkChainNode: [REPORT] [106254|PwBaseWorkChain|results]: work chain completed after 1 iterations
    05/09/2025 09:42:47 AM <99094> aiida.orm.nodes.process.workflow.workchain.WorkChainNode: [REPORT] [106254|PwBaseWorkChain|on_terminated]: remote folders will not be cleaned
