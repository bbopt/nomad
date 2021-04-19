.. _installation:

Installation
============

On Linux and Mac OS X, NOMAD can be compiled using *CMake*, a tool to manage building of source code.

The minimum version of *CMake* is 3.14. Older versions should trigger an error.

A recent C++ compiler is also required. *CMake* will detect which compiler is available.

.. warning:: Some older version of *CMake* do not trigger an error. If the ``cmake`` commands fail, check the version manually on the command line

  ::

    cmake --version

  **The minimum acceptable version is 3.14.**



.. note:: If the version of *CMake* is older than 3.14 or if you do not have *CMake* installed, you
   can follow the procedure given at `cmake.org <https://cmake.org/install/>`_.

   Alternatively, for Mac OSX, *CMake* can be installed on the command line using package manager `MacPorts <https://www.macports.org/>`_ or `Homebrew <http://brew.sh/>`_.


The procedure has the three following steps: **configuration, building and installation**.


.. _cmake_configuration:

Configuration using provided ``CMakeLists.txt`` files
"""""""""""""""""""""""""""""""""""""""""""""""""""""

On the command line in the ``$NOMAD_HOME`` directory::

  cmake -S . -B build/release

This command creates the ``CMake files`` and directories for building (``-B``) in ``build/release``. The source (``-S``) ``CMakeLists.txt`` file is in the ``$NOMAD_HOME`` directory.

To enable time stats build::

  cmake -DTIME_STATS=ON -S . -B build/release


To enable interfaces (C and Python) building::

  cmake -DBUILD_INTERFACES=ON -S . -B build/release


*OpenMP* is used for parallelization of evaluations. *CMake* will detect if *OpenMP* is available by default. To forcefully deactivate compilation with OpenMP::

  cmake -DTEST_OPENMP=OFF -S . -B build/release

Build
"""""

Build the libraries and applications::

  cmake --build build/release

Option ``--parallel xx`` can be added for faster build

Install
"""""""

Copy binaries and headers in build/release/[bin, include, lib] and in the examples/tests directories::

  cmake --install build/release

  or

  cmake --install build/debug

The executable ``nomad`` will installed into the directory::

  build/release/bin/

  or

  build/debug/bin/


Bulding for debug version
"""""""""""""""""""""""""

The procedure to build the ``debug`` version is the following. On the command line in the ``$NOMAD_HOME`` directory::

  cmake -S . -B build/debug -D CMAKE_BUILD_TYPE=Debug

  cmake --build build/debug

  cmake --install build/debug


Use another compiler
""""""""""""""""""""

The environment variables ``CC`` and ``CXX`` can be used to select the compiler.

.. note:: ``Clang`` is the default compiler for Mac OSX using XCode. Users of MAC OSX can install ``GCC`` compilers using `MacPorts <https://www.macports.org/>`_ or `Homebrew <http://brew.sh/>`_.


Testing installation
^^^^^^^^^^^^^^^^^^^^

Once building and installation have been performed some tests can be performed.
By default the examples are built and can be tested::

  cd build/release
  ctest

Option ``--parallel xx`` can be added for faster execution.
The log of the tests can be found in ``$NOMAD_HOME/build/release/Testing/Temporary``.
