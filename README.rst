==========================
Roodmus plugin for Scipion
==========================

Plugin to run `Roodmus <https://github.com/ccpem/roodmus>`_ software in Scipion.

All the information regarding the software is avalaible in the following `paper <https://www.biorxiv.org/content/10.1101/2024.04.29.590932v1>`_.

Roodmus is a benchmarking tool that allows for the analysis of cryo-EM reconstruction methods using ground truth data. Roodmus uses a molecular dynamics (MD) trajectory as input and samples this trajectory to generate a data set with encoded heterogeneity as defined by the MD. Synthetic micrographs are generated in the package using the `Parakeet <https://github.com/rosalindfranklininstitute/parakeet>`_ simulation software.

==========================
Plugin installation
==========================

There are two different ways to install Roodmus plugin:

Production mode
__________________________

.. code-block::

    scipion3 installp -p scipion-em-roodmus

Development mode
__________________________

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-roodmus.git
    scipion3 installp -p scipion-em-roodmus --devel
