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

Handling crashes in library mode requires special attention to isolate the part of code that may generate crashes. And, in general, using the library mode will require more computer programming than the batch mode. However, the library mode offers more options and flexibility for blackbox integration and management of optimization (see :ref:`library_mode`).

The different steps for solving your problem in batch mode are:

*  Create a directory for your problem. The problem directory is where the NOMAD command is executed. It is a convenient place to put the blackbox executable, the parameters file and the output files, but those locations can be customized.

*  Create your blackbox evaluation, which corresponds  to a program (a binary executable or a script). This program can be located in the problem directory or not.  This program outputs the objectives and the constraints for a given design vector. If you already have a blackbox program in a certain format, you need to interface it with a wrapper program to match  NOMAD specifications (see :ref:`getting_started` for blackbox basics).

*  Create a parameters file, for example ``param.txt``. This file can be located in the problem directory or not (see :ref:`basic_parameters_description` for more details).

* In the problem directory, start the optimization with a command like::

  $NOMAD_HOME/bin/nomad param.txt



.. _basic_parameters_description:

Basic parameters description
----------------------------

This section describes the basic parameters for the optimization problem definition, the algorithmic parameters and the parameters to manage output information. Additional information can be obtained by executing the command::

  $NOMAD_HOME/bin/nomad -h

to see all parameters, or::

  $NOMAD_HOME/bin/nomad -h PARAM_NAME

for a particular parameter.

The remaining content of a line is ignored after the character ``#``. Except for the file names, all strings and parameter names are case insensitive: ``DIMENSION 2`` is the same as ``Dimension 2``. File names refer to files in the problem directory. To indicate a file name containing spaces, use quotes ``"name"`` or ``'name'``. These names may include directory information relatively to the problem directory. The problem directory will be added to the names, unless the ``$`` character is used in front of the names. For example, if a blackbox executable is run by the command ``python script.py``, define parameter ``BB_EXE "$python script.py"``.

Some parameters consists of a list of variable indices taken from 0 to :math:`n-1`  (where :math:`n` is the number of variables). Variable indices may be entered individually or as a range  with format ``i-j``. Character ``*`` may be used to replace 0 to :math:`n-1`. Other parameters require arguments of type boolean: these values  may be entered with the strings ``yes``, ``no``, ``y``, ``n``,  ``0``, or ``1``. Finally, some parameters need vectors as arguments,  use ``(v1 v2 ... vn)`` for those. The strings ``-``, ``inf``, ``-inf`` or ``+inf``  are accepted to enter undefined real values  (NOMAD considers :math:`\pm \infty` as an undefined value).

Parameters are classified into problem, algorithmic and output parameters, and provided in what follows. The advanced functionalities of NOMAD are presented in :ref:`advanced_functionalities`.

.. _problem_parameters:

Problem parameters
^^^^^^^^^^^^^^^^^^

.. csv-table:: Basic problem parameters
   :header: "Name", "Argument", "Short description", "Default"
   :widths: 20,7,20,7

   :ref:`BB_EXE <bb_exe>`, list of strings, blackbox executables (required in batch mode) , Empty string
   :ref:`BB_INPUT_TYPE <bb_input_type>`, list of types, blackbox input types ,  ``* R`` (all real)
   :ref:`BB_OUTPUT_TYPE <bb_output_type>`, list of types , blackbox output types (required) , ``OBJ``
   DIMENSION, integer, :math:`n` the number of variables (required), 0
   :ref:`LOWER_BOUND <bounds>`, array of doubles  , lower bounds , none
   :ref:`UPPER_BOUND <bounds>`, array of doubles, upper bounds, none


.. _bb_exe:

``BB_EXE``
""""""""""

In batch mode, ``BB_EXE`` indicates the names of the blackbox executables.

A single string may be given if a single blackbox is used and gives several outputs. It is also possible to indicate several blackbox executables.

A blackbox program can return more than one function :ref:`BB_OUTPUT_TYPE <bb_output_type>`::

  BB_EXE         bb.exe             # defines that `bb.exe' is an
  BB_OUTPUT_TYPE OBJ EB EB          # executable with 3 outputs


A mapping between the names of the blackbox programs and the ``BB_OUTPUT_TYPE`` may be established to identify which function is returned by which blackbox::

  BB_EXE         bb1.exe bb2.exe    # defines two blackboxes
  BB_OUTPUT_TYPE OBJ     EB         # `bb1.exe' and `bb2.exe'
                                    # with one output each

Blackbox program names can be repeated to establish more complex mapping::

  BB_EXE   bb1.exe bb2.exe bb2.exe  # defines TWO blackboxes
                                    # NO duplication if names are repeated
  BB_OUTPUT_TYPE EB OBJ PB          # bb1.exe  has one output
                                    # bb2.exe  has two outputs
                                    # bb1.exe is executed first.
                                    #!! If EB constraint is feasible then
                                    #!!        bb2.exe is executed.
                                    #!! If EB constraint not feasible then
                                    #!!      bb2.exe is not launched.


A path can precede the blackbox program but spaces are not accepted in the path::

  BB_EXE "dir_of_blackbox/bb.exe"

To prevent NOMAD from adding a path,  the special character ``$``  should be put in front of a command::

  BB_EXE "$python bb.py"          # the blackbox is a python
                                  # script: it is run with
                                  # command
                                  # `python PROBLEM_DIR/bb.py'

Or::

  BB_EXE "$nice bb.exe"           # to run bb.exe
                                  # in nice mode on X systems



.. _bb_input_type:

``BB_INPUT_TYPE``
"""""""""""""""""

This parameter indicates the types of each variable.
It may be defined once with a list of :math:`n` input types with format  ``(t1 t2 ... tn)`` or `` * t``.
Input types ``t`` are values in ``R, B, I``. ``R`` is for real/continuous variables, ``B`` for binary variables, and ``I`` for integer variables.
The default type is ``R``. See also :doc:`ListOfParameters`.

.. note:: Categorical variables are not yet supported in NOMAD 4 but are available in NOMAD 3.




.. _bb_output_type:

``BB_OUTPUT_TYPE``
""""""""""""""""""

This parameter characterizes the values supplied by the blackbox, and in particular tells how constraint values are to be treated. The arguments are a list of :math:`m` types, where :math:`m` is the number of outputs of the  blackbox. At least one of these values must correspond to the objective function that NOMAD minimizes. Currently, NOMAD 4 only supports single objective problem (NOMAD 3 can handle bi-objective). Other values typically are constraints of the form :math:`c_j(x) \leq 0`, and the blackbox  must display the left-hand side of the constraint with this format.

.. note:: A terminology is used to describe the different types of constraints [AuDe09a]_

  * ``EB`` constraints correspond to constraints that need to be always satisfied (*unrelaxable constraints*). The technique used to deal with those is the **Extreme Barrier** approach, consisting in simply rejecting the  infeasible points.

  * ``PB`` and ``F`` constraints correspond to constraints that need to be satisfied only at the solution, and not necessarily at intermediate points (*relaxable constraints*). More precisely, ``F`` constraints are treated with the **Filter** approach [AuDe04a]_,  and ``PB`` constraints are treated with the **Progressive Barrier**  approach [AuDe09a]_.

  * There may be another type of constraints, the *hidden constraints*, but these only  appear inside the blackbox during an execution, and thus they   cannot be indicated in advance to NOMAD  (when such a constraint is violated, the evaluation simply fails and the point  is not considered).

  If the user is not sure about the nature of its constraints, we suggest using the keyword ``CSTR``, which corresponds by default to ``PB`` constraints.

All the types are:

+---------------+-------------------------------------------------------+
| ``CNT_EVAL``  |  Must be 0 or 1: count or not the blackbox evaluation |
+---------------+-------------------------------------------------------+
| ``EB``        | Constraint treated with **Extreme Barrier**           |
|               | (infeasible points are ignored)                       |
+---------------+-------------------------------------------------------+
| ``F``         | Constraint treated with **Filter** approach           |
+---------------+-------------------------------------------------------+
| ``NOTHING``   | The output is ignored                                 |
| ``EXTRA_O``   |                                                       |
| ``-``         |                                                       |
+---------------+-------------------------------------------------------+
| ``OBJ``       | Objective value to be minimized                       |
+---------------+-------------------------------------------------------+
| ``PB``        | Constraint treated with **Progressive Barrier**       |
| ``CSTR``      |                                                       |
+---------------+-------------------------------------------------------+


Please note that ``F`` constraints are not compatible with ``CSTR`` or ``PB``. However, ``EB`` can be combined  with ``F``, ``CSTR`` or ``PB``.





.. _bounds:

``LOWER_BOUND`` and ``UPPER_BOUND``
"""""""""""""""""""""""""""""""""""

.. warning:: NOMAD is 0 based :math:`\rightarrow` The first variable has a 0 index.

Parameters ``LOWER_BOUND`` and ``UPPER_BOUND`` are used to define bounds on variables. For example, with :math:`n=7`::

  LOWER_BOUND  0-2  -5.0
  LOWER_BOUND  3     0.0
  LOWER_BOUND  5-6  -4.0
  UPPER_BOUND  0-5   8.0


is equivalent to::

  LOWER_BOUND ( -5 -5 -5 0 - -4 -4 ) # `-' or `-inf' means that x_4
                                     # has no lower bound
  UPPER_BOUND (  8 8 8 8 8 8 inf )   # `-' or `inf' or `+inf' means
                                     # that x_6 has no upper bound.

Each of these two sequences define the following bounds

.. math::

  -5 ~ \leq x_0 \leq ~ 8 \\
  -5 ~ \leq x_1 \leq ~ 8 \\
  -5 ~ \leq x_2 \leq ~ 8 \\
   0 ~ \leq x_3 \leq ~ 8 \\
            x_4 \leq ~ 8 \\
  -4 ~ \leq x_5 \leq ~ 8 \\
  -4 ~ \leq x_6 \qquad \\


.. _algorithmic_parameters:

Algorithmic parameters
^^^^^^^^^^^^^^^^^^^^^^

.. csv-table:: Basic algorithmic parameters
   :header: "Name", "Argument", "Short description", "Default"
   :widths: 20,7,20,7

   :ref:`DIRECTION_TYPE <direction_type>` , direction type, type of directions for the poll , ``ORTHO 2N``
   F_TARGET , real :math:`t` , NOMAD terminates if :math:`f(x_k) \leq t` for the objective function, none
   :ref:`INITIAL_MESH_SIZE <initial_mesh_size>` , array of doubles, :math:`\delta_0` [AuDe2006]_ , none
   :ref:`INITIAL_FRAME_SIZE <initial_mesh_size>` , array of doubles, :math:`\Delta_0` [AuDe2006]_ , ``r0.1`` or based on ``X0``
   LH_SEARCH , 2 integers: ``p0`` and ``pi`` , **LH (Latin-Hypercube) search** (``p0``: initial and ``pi``: iterative) , none
   MAX_BB_EVAL , integer, maximum number of blackbox evaluations, none
   MAX_TIME , integer , maximum wall-clock time (in seconds) , none
   :ref:`TMP_DIR <tmp_dir>` , string, temporary directory for blackbox i/o files , problem directory
   :ref:`X0 <x0>` , point , starting point(s) , best point from a cache file or from an initial **LH search**


.. _direction_type:

``DIRECTION_TYPE``
""""""""""""""""""

This parameter defines the type of directions for *Mads* *Poll* step. The possible arguments are:

.. csv-table:: Direction types
   :widths: 6,20

   ``ORTHO 2N``, "OrthoMADS, 2n. This corresponds to the original *Ortho Mads* algorithm [AbAuDeLe09]_ with :math:`2n` directions."
   ``ORTHO N+1 NEG``, "OrthoMADS, n+1, with ((n+1)th dir = negative sum of the first n dirs) [AuIaLeDTr2014]_"
   ``N+1 UNI``, "MADS with n+1, using :math:`n+1` uniformly distributed directions."
   ``SINGLE``, "A single direction is produced"
   ``DOUBLE``, "Two opposite directions are produced"

Multiple direction types may be chosen by specifying ``DIRECTION_TYPE`` several times.


.. _initial_mesh_size:

``INITIAL_MESH_SIZE`` and ``INITIAL_FRAME_SIZE``
""""""""""""""""""""""""""""""""""""""""""""""""

The *Poll* step initial frame size :math:`\Delta_0` is decided by ``INITIAL_FRAME_SIZE``. In order to achieve the scaling between variables, NOMAD considers the frame size parameter for each variable independently. The initial mesh size parameter ``\delta_0`` is decided based on ``\Delta_0``. ``INITIAL_FRAME_SIZE`` may be entered with the following formats::

  INITIAL_FRAME_SIZE     d0 (same initial mesh size for all variables)
  INITIAL_FRAME_SIZE     (d0 d1 ... dn-1) (for all variables ``-`` may be used,  and defaults will be considered)
  INITIAL_FRAME_SIZE i   d0 (initial mesh size provided for variable ``i`` only)
  INITIAL_FRAME_SIZE i-j d0 (initial mesh size provided for variables ``i`` to ``j``)

The same logic and format apply for providing the ``INITIAL_MESH_SIZE``, ``MIN_MESH_SIZE`` and ``MIN_FRAME_SIZE``.


.. _tmp_dir:

``TMP_DIR``
"""""""""""

If NOMAD is installed on a network file system, with the batch mode use,  the cost of read/write files  will be high if no local temporary directory is defined.  On \linux/\unix/\osx\ systems, the directory ``/tmp`` is  local and we advise the user to define ``TMP_DIR /tmp``.


.. _x0:

``X0``
""""""

Parameter ``X0`` indicates the starting point of the algorithm. Several starting points may be proposed by entering this parameter several times. If no starting point is indicated, NOMAD considers the best evaluated point from an existing cache file (parameter ``CACHE_FILE``)  or from an initial *Latin-Hypercube search* (argument ``p0`` of ``LH_SEARCH``).

The ``X0`` parameter may take several types of arguments:

* A string indicating an existing cache file, containing several points (they can be already evaluated or not). This file may be the same as the one indicated with ``CACHE_FILE``. If so, this file will be updated during the program execution, otherwise the file will not be modified.

* A string indicating a text file containing the coordinates of one or several points (values are separated by spaces or line breaks).

* :math:`n` real values with format ``(v0 v1 ... vn-1)``.

* ``X0`` keyword plus integer(s) and one real

::

  X0 i v: (i+1)th coordinate set to v.
  X0 i-j v: coordinates i to j set to v.
  X0 * v: all coordinates set to v.

* One integer, another integer (or index range) and one real: the same as above except that the first integer k refers to the (k+1)th starting point.

The following example with :math:`n=3` corresponds to the two starting points :math:`( 5~0~0)` and :math:`(-5~1~1)`::

  X0 * 0.0
  X0 0 5.0
  X0 1 * 1.0
  X0 1 0 -5.0

.. _output_parameters:

Output parameters
^^^^^^^^^^^^^^^^^

.. csv-table:: Basic output parameters
   :header: "Name", "Argument", "Short description", "Default"
   :widths: 20,7,20,7

   CACHE_FILE , string ,  "cache file; if the file does not exist, it will be created", none
   DISPLAY_ALL_EVAL , bool , if ``yes`` all points are displayed with ``DISPLAY_STATS`` and ``STATS_FILE`` , ``no``
   DISPLAY_DEGREE , integer in [0 ; 3] or a string with four digits (see online help) , 0 no display and 3 full display ,  ``2``
   :ref:`DISPLAY_STATS <display_stats>` , list of strings ,  what information is displayed at each success , ``BBE OBJ``
   HISTORY_FILE , string ,  file containing all trial points with format ``x1 x2 ... xn bbo1 bbo2 ... bbom`` on each line ,  none
   SOLUTION_FILE, string , file to save the best feasible incumbent point in a simple format (SOL BBO) , none
   :ref:`STATS_FILE <display_stats>`  , string ``file_name`` + list of strings , the same as ``DISPLAY_STATS`` but for a display into file ,  none


.. _display_degree:

``DISPLAY_DEGREE``
""""""""""""""""""

Four different levels of display can be set via the parameter ``DISPLAY_DEGREE``. The ``DISPLAY_MAX_STEP_LEVEL`` can be used to control the number of steps displayed. To control the display of the **Models**, a ``QUAD_MODEL_DISPLAY`` and a ``SGTELIB_MODEL_DISPLAY`` are available. More information on these parameters can be obtained with online documentation: ``$NOMAD_HOME/bin/nomad -h display``


.. _display_stats:

``DISPLAY_STATS`` and ``STATS_FILE``
""""""""""""""""""""""""""""""""""""

These parameters display information each time a new feasible incumbent (i.e. a new best solution) is found.  ``DISPLAY_STATS`` is used to display at the standard output and ``STATS_FILE`` is used to write into a file.  These parameters need a list of strings as argument, **without any quotes**.  These strings may include the following keywords:

+----------+-------------------------------------+
| ``BBE``  |  The number of blackbox evaluations |
+----------+-------------------------------------+
| ``BBO``  |  The blackbox outputs               |
+----------+-------------------------------------+
| ``OBJ``  |  The objective function value       |
+----------+-------------------------------------+
| ``SOL``  |  The current feasible iterate       |
+----------+-------------------------------------+

.. note::
  More display options are available. Check the online help: ``$NOMAD_HOME/bin/nomad -h display_stats``


.. topic:: References

  .. [AuDe04a] C. Audet and J.E. Dennis, Jr. A pattern search filter method for nonlinear programming without derivatives. *SIAM Journal on Optimization*, 14(4):980–1010, 2004.
  .. [AuIaLeDTr2014] C. Audet and A. Ianni and S. Le~Digabel and C. Tribes. Reducing the Number of Function Evaluations in Mesh Adaptive Direct Search Algorithms. *SIAM Journal on Optimization*, 24(2):621-642, 2014.
