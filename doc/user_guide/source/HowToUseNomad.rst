.. _basic_nomad_usage:

NOMAD usage
===========

This chapter describes how to use NOMAD for solving blackbox optimization problems. Functionalities of NOMAD that are considered more advanced such as parallel evaluations are presented in :ref:`advanced_functionalities`.

.. note::

  New users are encouraged to first read :ref:`getting_started` to understand the basics of NOMAD utilization.


.. note::
  Many examples are provided in ``$NOMAD_HOME/examples`` with typical optimization outputs.

Batch mode is presented first, followed by a description of the basic parameters to setup and solve the majority of optimization problems that NOMAD can handle. The library mode is described in :ref:`library_mode`.

NOMAD should be cited with references [AuCo04a]_ and [AuLeRoTr2021]_. Other relevant papers by the developers are accessible through  the NOMAD website  `<http://www.gerad.ca/nomad>`_.

.. topic:: References

  .. [AuCo04a] M.A. Abramson, C. Audet, G. Couture, J.E. Dennis, Jr., S. Le Digabel, V. Rochon Montplaisir, and C. Tribes. The NOMAD project. Software available at `<https://www.gerad.ca/nomad>`_, 2021.

Optimization in batch mode
--------------------------

The batch mode allows to separate the evaluation of the objectives and constraints by the blackbox program from NOMAD executable. This mode has the advantage that if your blackbox program crashes, it will not affect NOMAD: The point that caused this crash will simply be tagged as a blackbox failure.


.. _basic_parameters_description:

Basic parameters description
----------------------------





.. _library_mode:

Optimization in library mode
----------------------------

...

Python interface
----------------

...
