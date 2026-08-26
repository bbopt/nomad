![NOMAD logo](https://www.gerad.ca/system/assets/000/001/537/1537.Nomad_original.jpg)

# NOMAD

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](LICENSE)

## Table of Contents

- [Description](#description)
- [Web Pages](#web-pages)
- [Contact](#contact)
- [Version Warning](#version-warning)
- [How to Cite](#how-to-cite)
- [Compilation (Release)](#compilation-release)
- [Compilation (Debug)](#compilation-debug)
- [Examples of Optimization](#examples-of-optimization)
- [Authors and License](#authors-and-license)
- [Contributing](#contributing)

## Description

NOMAD is a C++ implementation of the Mesh Adaptive Direct Search (MADS)
algorithm, designed for constrained optimization of black-box functions.

The algorithms implemented are based on the book
"Derivative-Free and Blackbox Optimization", by Charles Audet and Warren Hare,
Springer 2017. A second edition of the book has been released in 2026.

## Web Pages

- <https://www.gerad.ca/nomad/>
- <https://github.com/bbopt/nomad>

## Contact

- Polytechnique Montreal - GERAD
- C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada
- e-mail: nomad@gerad.ca

## Version Warning

This repository is for NOMAD 4. The previous version, NOMAD 3,
is not on GitHub. NOMAD 3 will be deprecated in the future.

NOMAD 4 is similar in usage to NOMAD 3.

NOMAD 4 has a new software architecture, uses OpenMP to run
evaluations in parallel, and also has some new functionalities.

## How to Cite

Please cite NOMAD 4 with reference:

> C. Audet, S. Le Digabel, V. Rochon Montplaisir, and C. Tribes.
> Algorithm 1027: NOMAD version 4: Nonlinear optimization with the
> MADS algorithm. ACM Transactions on Mathematical Software
> Volume 48, Issue 3, Article No.: 35, pp 1–22
> <https://doi.org/10.1145/3544489>

## Compilation (Release)

On Linux, Unix, Windows and Mac OS X, NOMAD can be compiled using CMake.
The minimum version of CMake is 3.14. Older versions will trigger
an error. A recent C++ compiler is also required.

Before starting the procedure, we recommend setting the environment variable
`$NOMAD_HOME` with the path where NOMAD has been copied. For Linux and OSX,

```bash
export NOMAD_HOME=/home/myUserName/PathToNomad
```

For Windows, add an environment variable `%NOMAD_HOME%` containing the path.
The remaining text uses the `$NOMAD_HOME` environment variable.

The procedure is the following. On the command line in the
`$NOMAD_HOME` directory:

```bash
cmake -S . -B build/release
```
---> Create the CMake files and directories for building (`-B`) in `build/release`.
     The source (`-S`) `CMakeLists.txt` file is in the `$NOMAD_HOME` directory.

- To deactivate compilation with OpenMP:
  ```bash
  cmake -DTEST_OPENMP=OFF -S . -B build/release
  ```

- To enable *C interface* building:
  ```bash
  cmake -DBUILD_INTERFACE_C=ON -S . -B build/release
  ```

- To enable *Matlab* interface building:
  ```bash
  cmake -DBUILD_INTERFACE_MATLAB=ON  -DTEST_OPENMP=OFF -S . -B build/release
  ```

  - ! Before proceeding, have a look into
    `$NOMAD_HOME/interfaces/Matlab_MEX/readme.txt`

  - ! Extra flags might be required to prevent CMake errors.

  - ! Compiler version and Matlab version need to be compatible;
    Check <https://www.mathworks.com/support/requirements/supported-compilers.html>

  - ! The Matlab interface will not be built if OpenMP is enabled.

  - ! An extra `addpath` Matlab command must be done to have access
    to nomad Mex binaries:
    ```matlab
    addpath(strcat(getenv('NOMAD_HOME'),'/build/release/bin'))
    ```

- To enable *Python* interface (PyNomad) building:
  ```bash
  cmake -DBUILD_INTERFACE_PYTHON=ON -S . -B build/release
  ```

  - ! Before proceeding, have a look into
    `$NOMAD_HOME/interfaces/PyNomad/readme.txt`

  - ! More details are provided in `$NOMAD_HOME/interfaces/PyNomad/readme.txt`

  - ! Building requires Cython. Cython can be obtained with
    Anaconda distribution platform.

  - ! On *Windows*, using Visual Studio, see the user guide or the README
    to properly manage x86/x64 building of binaries.

- To enable *Java* interface building (with Swig):
  ```bash
  cmake -DBUILD_INTERFACE_JAVA=ON -S . -B build/release
  ```

For *Windows*:
```bash
cmake --build build/release --config Release
```
For *OSX* and *Linux*:
```bash
cmake --build build/release
```
---> Build all the libraries and applications

- Option `--parallel xx` can be added for faster build.
- The option `--config Release` should be used on *Windows*
  multi-configuration build environment (VisualStudio) to build only
  Release configuration. The default configuration is Debug.
  The same option should be used for *OSX* when using a *Xcode* project.

For *Windows*:
```bash
cmake --install build/release --config Release
```
For *OSX* and *Linux*:
```bash
cmake --install build/release
```
---> Copy binaries and headers in `build/release/[bin, include, lib]`
     and in the examples/tests directories.

By default, the executable "nomad" will be installed into the directory:
`build/release/bin/` (`build/debug/bin/` when in debug mode). A symbolic link
is added in the bin directory (not functional for windows).

It is possible to build only a single application in its working directory
(with `NOMAD_HOME` environment variable properly set):

```bash
cd $NOMAD_HOME/examples/basic/library/example1
cmake --build $NOMAD_HOME/build/release --target example1_lib.exe
cmake --install $NOMAD_HOME/build/release
```

## Compilation (Debug)

The procedure to build the debug version is the following.
On the command line in the `$NOMAD_HOME` directory:

```bash
cmake -S . -B build/debug -D CMAKE_BUILD_TYPE=Debug
```
---> On *Windows*, all 4 configurations are configured
     (Debug, RelWithDebugInfo, MinSizeRel, Release); flag
     `CMAKE_BUILD_TYPE` is ignored.

For *Windows*:
```bash
cmake --build build/debug --config Debug
```
For *OSX* and *Linux*:
```bash
cmake --build build/debug
```
---> Build the libraries and applications

- Option `--parallel xx` can be added for faster build.

For *Windows*:
```bash
cmake --install build/debug --config Debug
```
For *OSX* and *Linux*:
```bash
cmake --install build/debug
```
---> Copy binaries and headers in `build/debug/[bin, include, lib]`
     and in the examples/tests directories

## Examples of Optimization

### Batch Mode

There are examples in batch mode in `$NOMAD_HOME/examples/basic/batch/`.
In each directory, the blackbox functions (usually named bb) are compiled
by default. The problem may be resolved using NOMAD and the parameter file:

```bash
$NOMAD_HOME/build/release/bin/nomad param.txt
```

For convenience, the path to `$NOMAD_HOME/build/release/bin` directory
can be added to the `$PATH` environment variable. For *Windows*, this is
achieved by setting the parameters for environment variable `%PATH%` to
`%NOMAD_HOME\build\release\bin\`

### Library Mode

There are examples in library mode in `$NOMAD_HOME/examples/basic/library/`.
In each directory, the executable may be compiled when building
the NOMAD application. The problems may be resolved by execution,
for instance:

```bash
./example_lib.exe
```

**IMPORTANT**: Library mode examples with *Windows* require to set
the `%PATH%` environment variable (see above). Otherwise, the executables
cannot find NOMAD dlls.

## Authors and License

NOMAD - Version 4 has been created and developed by Viviane Rochon Montplaisir
and Christophe Tribes (Polytechnique Montreal).

The copyright of NOMAD - version 4 is owned by Charles Audet, Sebastien Le Digabel,
Viviane Rochon Montplaisir and Christophe Tribes (Polytechnique Montreal).

NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada, NSERC
(Natural Sciences and Engineering Research Council of Canada), InnovÉÉ
(Innovation en Énergie Électrique) and IVADO (The Institute for Data Valorization).

NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,
Christophe Tribes and Viviane Rochon Montplaisir, and was funded by AFOSR and
Exxon Mobil.

NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,
Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and Exxon Mobil.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version. See the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
