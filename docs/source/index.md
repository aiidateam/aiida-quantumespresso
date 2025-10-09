---
myst:
  substitutions:
    README.md of the repository: '`README.md` of the repository'
    aiida-core documentation: '`aiida-core` documentation'
    aiida-quantumespresso: '`aiida-quantumespresso`'
---

```{toctree}
:hidden: true
:caption: Getting Started

get_started/installation
get_started/quick_start
tutorials/index
```

```{toctree}
:hidden: true
:caption: How to

howto/customize_inputs
howto/calculations/index
howto/workflows/index
```

```{toctree}
:hidden: true
:caption: Topic guides
topics/calculations/index
topics/workflows/index
```

```{toctree}
:hidden: true
:caption: Reference
reference/cli/index
```

```{toctree}
:hidden: true

developer
```

# AiiDA Quantum ESPRESSO

An AiiDA plugin package to integrate the [Quantum ESPRESSO](http://www.quantumespresso.org) software suite.
Compute a variety of material properties with the popular open source DFT code with automatic data provenance provided by AiiDA.
Geometry optimizations, ground-state electronic structure, band structures, phonons, and much more.

[![PyPI version](https://badge.fury.io/py/aiida-quantumespresso.svg)](https://badge.fury.io/py/aiida-quantumespresso)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/aiida-quantumespresso.svg)](https://pypi.python.org/pypi/aiida-quantumespresso)
[![Build Status](https://github.com/aiidateam/aiida-quantumespresso/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/aiidateam/aiida-quantumespresso/actions)
[![Docs status](https://readthedocs.org/projects/aiida-quantumespresso/badge)](http://aiida-quantumespresso.readthedocs.io/)

______________________________________________________________________


::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fa}`info-circle;mr-1` Installation
:text-align: center
:shadow: md

Instructions to install, configure and setup the plugin package.

+++

```{button-ref} get_started/installation
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:

To the installation guides
```
:::

:::{grid-item-card} {fa}`rocket;mr-1` Get started
:text-align: center
:shadow: md

Easy examples to take the first steps with the plugin package.

+++

```{button-ref} get_started/quick_start
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:

To the tutorial
```
:::
::::


## How to cite

If you use this plugin for your research, please cite the following work:

> Sebastiaan. P. Huber, Spyros Zoupanos, Martin Uhrin, Leopold Talirz, Leonid Kahle, Rico Häuselmann, Dominik Gresch, Tiziano Müller, Aliaksandr V. Yakutovich, Casper W. Andersen, Francisco F. Ramirez, Carl S. Adorf, Fernando Gargiulo, Snehal Kumbhar, Elsa Passaro, Conrad Johnston, Andrius Merkys, Andrea Cepellotti, Nicolas Mounet, Nicola Marzari, Boris Kozinsky, and Giovanni Pizzi, [*AiiDA 1.0, a scalable computational infrastructure for automated
    reproducible workflows and data provenance*](https://doi.org/10.1038/s41597-020-00638-4), Scientific Data **7**, 300 (2020)

> Martin Uhrin, Sebastiaan. P. Huber, Jusong Yu, Nicola Marzari, and Giovanni Pizzi, [*Workflows in AiiDA: Engineering a high-throughput, event-based
    engine for robust and modular computational workflows*](https://doi.org/10.1016/j.commatsci.2020.110086), Computational Materials Science **187**, 110086 (2021)

## Acknowledgements

We acknowledge support from:

:::{list-table}
:widths: 60 40
:class: logo-table
:header-rows: 0

* - The [NCCR MARVEL](http://nccr-marvel.ch/) funded by the Swiss National Science Foundation.
  - ![marvel](images/MARVEL.png)
* - The EU Centre of Excellence ["MaX – Materials Design at the Exascale"](http://www.max-centre.eu/) (Horizon 2020 EINFRA-5, Grant No. 676598).
  - ![max](images/MaX.png)
* - The [swissuniversities P-5 project "Materials Cloud"](https://www.materialscloud.org/swissuniversities)
  - ![swissuniversities](images/swissuniversities.png)

:::

[aiida]: http://aiida.net
[aiida quantum espresso tutorial]: https://aiida-tutorials.readthedocs.io/en/tutorial-qe-short/
[aiida-core documentation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html
[aiida-quantumespresso]: https://github.com/aiidateam/aiida-quantumespresso
[aiidalab demo cluster]: https://aiidalab-demo.materialscloud.org/
[quantum espresso]: http://www.quantumespresso.org
[quantum mobile]: https://quantum-mobile.readthedocs.io/en/latest/index.html
