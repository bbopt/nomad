Nomad 4 User Guide
==================

.. warning::
   NOMAD is a blackbox optimization software.

   This user guide is specific to NOMAD 4.

   NOMAD 3 is still available. It will be replaced by NOMAD 4 in the future.

   Get NOMAD 3 and 4 at `<https://www.gerad.ca/nomad/>`_.


A general presentation of NOMAD is given in :ref:`introduction`.

**New users** of NOMAD should refer to

* :ref:`installation`
* :ref:`getting_started`

**Using NOMAD**

* Starting from :ref:`basic_nomad_usage`, all users can find ways to tailor problem definition, algorithmic settings and software output.

* Refer to :ref:`advanced_functionalities` and :ref:`tricks_of_the_trade` for specific problem solving.

Please cite NOMAD 4 with reference:

.. [AuLeRoTr2021] C. Audet, S. Le Digabel, V. Rochon Montplaisir, and C. Tribes.
   NOMAD version 4: Nonlinear optimization with the MADS algorithm.
   *ACM Transactions on Mathematical Software*, Submitted.

A complete introduction to derivative-free and blackbox optimization can be found in the textbook:

.. [AuHa2017] C. Audet and W. Hare.
    Derivative-Free and Blackbox Optimization.
    *Springer Series in Operations Research and Financial Engineering.*
    Springer International Publishing, Berlin, 2017. |dfo_book|

.. |dfo_book| image:: ../figs/livre_DFO_AuHa2017.png
               :width: 30 pt

.. toctree::
   :maxdepth: 1
   :caption: Introduction:

   Introduction

.. toctree::
   :maxdepth: 1
   :caption: Installation and tests:

   Installation

.. toctree::
   :maxdepth: 1
   :caption: First NOMAD steps:

   GettingStarted
   HowToUseNomad
   TricksOfTheTrade
   AdvancedFunctionalities


.. toctree::
   :maxdepth: 1
   :caption: Additional information:

   ReleaseNotes
   Appendix

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
