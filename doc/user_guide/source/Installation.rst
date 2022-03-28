.. _installation:

Installation
============

On Linux, Windows and Mac OS X, NOMAD can be compiled using *CMake*, a tool to manage building of source code.

The minimum version of *CMake* is 3.14. Older versions should trigger an error. *CMake* will detect which compiler is available.

A recent C++ compiler supporting C++14 is also required. The compilation has been tested on Linux with gcc 9.3.0, 10.1.0 and 11.1.0. The compilation has been tested on OSX with gcc Homebrew 9.3.0 and 11.2.0, and also Apple clang version 11.0.3 and 13.0.0. The compilation has been tested on Windows 8 and Windows 10 Enterprise with Microsoft Visual Studio 2019 (cl.exe 19.29.300038.1) and Microsoft Visual Studio 2017.


.. warning:: Some older version of *CMake* do not trigger an explicit error on the version number.
   If the ``cmake`` commands fail, check the version manually on the command line

  ::

    cmake --version

  **The minimum acceptable version is 3.14.**



.. note:: If the version of *CMake* is older than 3.14 or if you do not have *CMake* installed,
   we recommend to install *CMake* using a **package manager**. The other option is to
   follow the procedure given at `cmake.org <https://cmake.org/install/>`_ to obtain binaries.

   For Mac OSX, *CMake* can be installed on the command line using package manager `MacPorts <https://www.macports.org/>`_ or `Homebrew <http://brew.sh/>`_.

   For Linux, several package managers exist to automate the procedure.

   For Windows, an installer tool is available at `cmake.org/download <https://cmake.org/download/>`_. Please note that all commands are performed in the Windows Command Prompt windows of Visual Studio.


The NOMAD installation procedure has the three following steps: **configuration, building and installation**.

.. warning:: Before starting the procedure we recommend to set the environment variable ``$NOMAD_HOME`` with the path where NOMAD has been copied. For Linux and OSX,

  ::

    export NOMAD_HOME=/home/myUserName/PathToNomad

For Windows, add an environment variable ``%NOMAD_HOME%`` containing the path.


  The remaining of the documentation uses the ``$NOMAD_HOME`` environment variable.





.. _cmake_configuration:

1- Configuration using provided ``CMakeLists.txt`` files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

On the command line, in the ``$NOMAD_HOME`` directory::

  cmake -S . -B build/release


.. sidebar:: Building options

     To enable time stats build::

        cmake -DTIME_STATS=ON -S . -B build/release

     To enable C interface building::

        cmake -DBUILD_INTERFACE_C=ON -S . -B build/release

     To enable Matlab interface building::

        cmake -DBUILD_INTERFACE_MATLAB=ON -S . -B build/release

     To enable Python interface (PyNomad) building::

        cmake -DBUILD_INTERFACE_PYTHON=ON -S . -B build/release

     To disable *OpenMP* compilation::

       cmake -DTEST_OPENMP=OFF -S . -B build/release


This command creates the files and directories for building (``-B``) in ``build/release``. The source (``-S``) ``CMakeLists.txt`` file is in the ``$NOMAD_HOME`` directory.

The command can be modified to enable/disable some options (see side bar).

*OpenMP* is used for parallelization of evaluations. *CMake* will detect if *OpenMP* is available by default. To forcefully deactivate compilation with *OpenMP*, see option in side bar.

.. warning:: Sometimes, additional options must be passed during configuration to prevent CMake errors for interfaces. See details in :ref:`library_mode`.



2- Build
""""""""

Build the libraries and applications (Linux/OSX)::

  cmake --build build/release

For Windows, the default configuration is Debug. To obtain the Release version::

  cmake --build build/release --config Release

Option ``--parallel xx`` can be added for faster build

It is possible to build only a single application in its working directory::

  cd $NOMAD_HOME/examples/basic/library/example1
  cmake --build $NOMAD_HOME/build/release --target example1_lib.exe

3- Install
""""""""""

Copy binaries and headers in build/release/[bin, include, lib] and in the examples/tests directories::

  cmake --install build/release

Option --config Release should be used on Windows to install Release configuration.

The executable ``nomad`` will installed into the directory::

  $NOMAD_HOME/build/release/bin/

Additionally a symbolic link to ``nomad`` binary is available::

  $NOMAD_HOME/bin



Bulding for debug version
"""""""""""""""""""""""""

The procedure to configure, build and install the ``debug`` version is the following (linux/OSX). On the command line in the ``$NOMAD_HOME`` directory::

  cmake -S . -B build/debug -D CMAKE_BUILD_TYPE=Debug

  cmake --build build/debug

  cmake --install build/debug

On Windows, all 4 configurations are always build Debug, RelWithDebugInfo, MinSizeRel, Release); the flag CMAKE_BUILD_TYPE can be ignored.

Use another compiler
""""""""""""""""""""

The environment variables ``CC`` and ``CXX`` can be used to select the ``C`` and ``C++`` compilers.

.. note:: ``Clang`` is the default compiler for Mac OSX using XCode. But, *OpenMP* (used for parallel evaluations)
   support is disabled in *Clang* that come with *Xcode*.
   Users of Mac OSX can install and use another compiler to enable *OpenMP* support.
   For example, ``GCC`` compilers can be obtained using `MacPorts <https://www.macports.org/>`_ or `Homebrew <http://brew.sh/>`_.


Testing installation
====================

Once building **and installation** have been performed some tests can be performed.

The NOMAD binary can be tested::

  $NOMAD_HOME/bin/nomad -v

This should return the version number on the command line.

By default the examples are built and can be tested::

  cd build/release
  ctest

For Windows, the configuration must be provided: ``ctest -C Release``.

Please note that the tests will take several minutes. Option ``--parallel xx`` can be added for faster execution.
The log of the tests can be found in ``$NOMAD_HOME/build/release/Testing/Temporary``.
