.. _library_mode:

Optimization in library mode
----------------------------

The library mode allows to tailor the evaluation of the objectives and constraints within a  specialized executable that contains NOMAD shared object library. For example, it is possible to link your own code with the NOMAD library (provided during installation) in a single executable that can define and run optimization for your problem. Contrary to the batch mode, this has the disadvantage that a crash within the executable (for example during the evaluation of a point) will end the optimization unless a special treatment of exception is provided by the user. But, as a counterpart, it offers more options and flexibility for blackbox integration and optimization management (display, pre- and post-processing, multiple optimizations, user search, etc.).

The library mode requires additional coding and compilation before conducting optimization. First, we will briefly review the compilation of source code to obtain NOMAD binaries (executable and shared object libraries) and how to use library.  Then, details on how to interface your own code are presented.

Compilation of the source code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOMAD source code provided during installation are located in ``$NOMAD_HOME/src``.  Examples are provided in ``$NOMAD_HOME/examples/basic/library`` and ``$NOMAD_HOME/examples/advanced/library``.

The compilation procedure uses provided ``CMake`` files along with the sources.

In what follows it is supposed that you have a write access to the source codes directory. If it is not the case, please consider making a copy in a more convenient location.

Using NOMAD libraries
^^^^^^^^^^^^^^^^^^^^^

Using the routines that are in the pre-compiled NOMAD shared object library (\unix / \linux / \osx) or dll (Windows) with a \cpp\ program requires building an executable. This is illustrated on the example located in the directory::

  $NOMAD_HOME/examples/basic/library/example1

For this example, just one ``C++`` source file is used, but there could be a lot more.


Python interface
----------------

...
