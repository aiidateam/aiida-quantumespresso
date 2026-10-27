(tutorials)=

# Tutorials

:::{important}
Before you get started, make sure that you have:

- installed the `aiida-quantumespresso` package ([see instructions](#installation))
- configured the `pw.x` code ([see instructions](#installation-setup-code))
- installed the SSSP pseudopotential family ([see instructions](#installation-setup-pseudopotentials))
:::

```{toctree}
:hidden: true
:maxdepth: 1

magnetism
hubbard
```

:::{card}
:class-header: panel-header-text
:class-footer: tutor-footer
:link: quick-start
:link-type: ref
:margin: 4

{fa}`fa-regular fa-rocket` Quick start: Running your first `pw.x` calculation
^^^
If you haven't already, go through the quick start tutorial.
Here you'll learn the basics on how to run the Quantum ESPRESSO `pw.x` code with AiiDA.
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


:::{card}
:class-header: panel-header-text
:class-footer: tutor-footer
:link: tutorials-hubbard
:link-type: ref
:margin: 4

{fa}`fa-solid fa-fire` Hubbard corrections
^^^
Learn how to define the Hubbard parameters along with your structure, and run a DFT+_U_+_V_ calculation using the `pw.x` binary.
+++
::::{list-table}
:class: footer-table
:widths: 50 50
* - {fa}`fa-sharp fa-regular fa-clock` 30 min
  - {{ aiida_logo }} [Beginner]{.aiida-green}
::::
:::
