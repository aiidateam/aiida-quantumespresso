(tutorials)=

# Tutorials

:::{important}
Before you get started, make sure that you have:

- installed the `aiida-quantumespresso` package ([see instructions](#installation-installation))
- configured the `pw.x` code ([see instructions](#installation-setup-code))
- installed the SSSP pseudopotential family ([see instructions](#installation-setup-pseudopotentials))
:::

```{toctree}
:hidden: true
:maxdepth: 1

first_pw
magnetism
```

:::{card}
:class-header: panel-header-text
:class-footer: tutor-footer
:link: tutorials-pw-through-cli
:link-type: ref
:margin: 4

{fa}`fa-regular fa-rocket` Running your first `pw.x` calculation
^^^
Learn how to run the Quantum ESPRESSO `pw.x` with AiiDA, both from the command line and using the Python API.
+++
::::{list-table}
:class: footer-table
:widths: 50 50
* - {fa}`fa-sharp fa-regular fa-clock` 30 min
  - {{ aiida_logo }} [Beginner]{.aiida-green}
::::
:::


:::{card}
:class-header: panel-header-text
:class-footer: tutor-footer
:link: tutorials-magnetic-configurations
:link-type: ref
:margin: 4

{fa}`fa-solid fa-arrow-down-up-across-line` Magnetic configurations
^^^
Learn how to assign magnetic configurations to your structure, and retrieve the final magnetic configuration from a `pw.x` calculation.
+++
::::{list-table}
:class: footer-table
:widths: 50 50
* - {fa}`fa-sharp fa-regular fa-clock` 20 min
  - {{ aiida_logo }} [Beginner]{.aiida-green}
::::
:::
