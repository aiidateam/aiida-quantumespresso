
.. grid::
   :reverse:
   :gutter: 2 3 3 3
   :margin: 1 2 1 2

   .. grid-item::
      :columns: 12 4 4 4

      .. image:: images/logo_aiida_quantumespresso.png
         :width: 200px
         :class: sd-m-auto

   .. grid-item::
      :columns: 12 8 8 8
      :child-align: justify
      :class: sd-fs-5

      .. rubric:: AiiDA Quantum ESPRESSO

      An AiiDA plugin package to integrate the `Quantum ESPRESSO`_ software suite.
      Compute a variety of material properties with the popular open source DFT code with automatic data provenance provided by AiiDA.
      Geometry optimizations, ground-state electronic structure, band structures, phonons, and much more.

      **aiida-quantumespresso version:** |release|


------------------------------

.. grid:: 1 2 2 2
    :gutter: 3

    .. grid-item-card:: :fa:`rocket;mr-1` Get started
        :text-align: center
        :shadow: md

        Instructions to install, configure and setup the plugin package.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: installation/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the installation guides

    .. grid-item-card:: :fa:`info-circle;mr-1` Tutorials
        :text-align: center
        :shadow: md

        Easy examples to take the first steps with the plugin package.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: tutorials/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the tutorials

    .. grid-item-card:: :fa:`question-circle;mr-1` How-to guides
        :text-align: center
        :shadow: md

        Hands-on guides to achieve specific goals.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: howto/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the how-to guides

    .. grid-item-card:: :fa:`bookmark;mr-1` Topic guides
        :text-align: center
        :shadow: md

        Detailed background information on various concepts.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: topics/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the topic guides

    .. grid-item-card:: :fa:`cogs;mr-1` Reference guides
        :text-align: center
        :shadow: md

        Detailed reference guides on the application programming and command line interfaces.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: reference/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the reference guides

.. toctree::
   :maxdepth: 2
   :hidden:

   installation/index
   tutorials/index
   howto/index
   topics/index
   reference/index


***********
How to cite
***********

If you use this plugin for your research, please cite the following work:

.. highlights:: Sebastiaan. P. Huber, Spyros Zoupanos, Martin Uhrin, Leopold Talirz, Leonid Kahle, Rico Häuselmann, Dominik Gresch, Tiziano Müller, Aliaksandr V. Yakutovich, Casper W. Andersen, Francisco F. Ramirez, Carl S. Adorf, Fernando Gargiulo, Snehal Kumbhar, Elsa Passaro, Conrad Johnston, Andrius Merkys, Andrea Cepellotti, Nicolas Mounet, Nicola Marzari, Boris Kozinsky, and Giovanni Pizzi, |AiiDA main paper|_, Scientific Data **7**, 300 (2020)

.. highlights:: Martin Uhrin, Sebastiaan. P. Huber, Jusong Yu, Nicola Marzari, and Giovanni Pizzi, |AiiDA engine paper|_, Computational Materials Science **187**, 110086 (2021)


****************
Acknowledgements
****************

We acknowledge support from:

.. list-table::
    :widths: 60 40
    :class: logo-table
    :header-rows: 0

    * - The `NCCR MARVEL`_ funded by the Swiss National Science Foundation.
      - |marvel|
    * - The EU Centre of Excellence "`MaX – Materials Design at the Exascale`_" (Horizon 2020 EINFRA-5, Grant No. 676598).
      - |max|
    * - The `swissuniversities P-5 project "Materials Cloud"`_.
      - |swissuniversities|

.. |marvel| image:: images/MARVEL.png
    :width: 100%

.. |max| image:: images/MaX.png
    :width: 100%

.. |swissuniversities| image:: images/swissuniversities.png
    :width: 100%

.. |aiida-core documentation| replace:: ``aiida-core`` documentation
.. _aiida-core documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html

.. |aiida-quantumespresso| replace:: ``aiida-quantumespresso``
.. _aiida-quantumespresso: https://github.com/aiidateam/aiida-quantumespresso

.. _AiiDA Quantum ESPRESSO tutorial: https://aiida-tutorials.readthedocs.io/en/tutorial-qe-short/

.. _AiiDA: http://aiida.net
.. _Quantum ESPRESSO: http://www.quantumespresso.org
.. _Quantum Mobile: https://quantum-mobile.readthedocs.io/en/latest/index.html
.. _AiiDAlab demo cluster: https://aiidalab-demo.materialscloud.org/

.. |AiiDA main paper| replace:: *AiiDA 1.0, a scalable computational infrastructure for automated reproducible workflows and data provenance*
.. _AiiDA main paper: https://doi.org/10.1038/s41597-020-00638-4

.. |AiiDA engine paper| replace:: *Workflows in AiiDA: Engineering a high-throughput, event-based engine for robust and modular computational workflows*
.. _AiiDA engine paper: https://doi.org/10.1016/j.commatsci.2020.110086

.. _NCCR MARVEL: http://nccr-marvel.ch/
.. _MaX – Materials Design at the Exascale: http://www.max-centre.eu/
.. _`swissuniversities P-5 project "Materials Cloud"`: https://www.materialscloud.org/swissuniversities

.. |README.md of the repository| replace:: ``README.md`` of the repository
.. _README.md of the repository: https://github.com/aiidateam/aiida-quantumespresso/blob/develop/README.md
