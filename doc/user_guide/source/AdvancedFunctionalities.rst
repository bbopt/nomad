.. _advanced_functionalities:

Advanced functionalities
========================

Advanced parameters
-------------------

Advanced parameters are intended to setup optimization problems, algorithmic and output parameters when specific needs are present. Only a few advanced parameters are presented below; all advanced parameters can be obtained with ``$NOMAD_HOME -h advanced``. Also a complete list of parameters and a short description is available in :ref:``appendix_parameters``.

.. _fixed_variable:

``FIXED_VARIABLE``
""""""""""""""""""

This parameter is used to fix some variables to a value. This value is optional if at least one starting point is defined. The parameter may be entered with several types of arguments:

* A vector of :math:`n values with format ``(v0 v1 ... vn-1)``. Character ```-`` is used for free variables.

* An index range if at least one starting point has been defined. ``FIXED_VARIABLE i-j``: variables ``i`` to ``j`` are fixed to their initial (``i-j`` may be replaced by ``i`` only). See :ref`x0` for practical examples of index ranges.

.. _seed:

``SEED``
""""""""

The directions that NOMAD explores during the *Poll* phase are dependent upon the seed. The seed is used to generate a pseudo-random direction on a unit n-dimensional sphere. The user can change the sequence of directions by setting ``SEED`` to a positive integer or ``-1``. If ``-1`` or ``DIFF`` is entered the seed is different for each run (PID is used).

Other aspects of NOMAD may depend on a pseudo-random sequence of numbers depending on selected options: *LH Search* and *PSD Mads*.

.. _eval_opportunistic:

``EVAL_OPPORTUNISTIC``
""""""""""""""""""""""

The opportunistic strategy consists in terminating the evaluations of a list of trial points at a given step of the algorithm as soon as an improved value is found.

This strategy is decided with the parameter ``EVAL_OPPORTUNISTIC`` and applies to both the *Poll* and *Search* steps. Search with NOMAD help ``$NOMAD_HOME/bin/nomad -h OPPORTUNISTIC`` for more options.

When evaluations are performed by blocks (see :ref:`bloc_evaluations`) the opportunistic strategy applies after evaluating a block of trial points.

.. _variable_group:

``VARIABLE_GROUP``
""""""""""""""""""

By default NOMAD creates one group that combines all continuous, integer, and binary variables.

In batch mode, the ``VARIABLE_GROUP`` parameter followed by variable indices is used to explicitly form a group of variables. Each group of variable generates its own polling directions. The parameter may be entered several times to define more than one group of variables. Variables in a group may be of different types.

.. _quad_model_search:

``QUAD_MODEL_SEARCH`` and ``SGTELIB_MODEL_SEARCH``
""""""""""""""""""""""""""""""""""""""""""""""""""

The *Search* phase of the *MADS* algorithm can use models of the objectives and constraints that are constructed dynamically from all the evaluations made. By default, a quadratic model is used to propose new points to be evaluated with the blackbox. To disable the use of quadratic models, the parameter ``QUAD_MODEL_SEARCH`` can be set to ``no``.

Models from the *SGTELIB* library can be used by setting ``SGTELIB_MODEL_SEARCH`` to ``yes``. Many parameters are available to control *SGTELIB* models: ``$NOMAD_HOME/bin/nomad -h SGTELIB``.

.. _granularity:

``GRANULARITY``
"""""""""""""""

The *MADS* algorithm handles granular variables, i.e. variables with a controlled number of decimals. For real numbers the granularity is 0. For integers and binary variables the granularity is automatically set to one.

The possible syntaxes to specify the granularity of the variables are as follows:

* :math:`n` real values with format ``GRANULARITY (v0 v1 ... vn-1)``.

* ``GRANULARITY i-j v``: coordinates  ``i` to  ``j`` set to ``v``.

* ``GRANULARITY * v``: all coordinates set to ``v``.




.. _block_evaluations:

Blackbox evaluation of a block of trial points
----------------------------------------------


.. _parallel_evaluations:

Parallel evaluations
--------------------

...

.. _psd_mads:

PSD Mads
--------

...
