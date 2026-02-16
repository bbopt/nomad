###################################################################################
#                                                                                 #
#                                    README                                       #
#                                                                                 #
#---------------------------------------------------------------------------------#
#  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                #
#                                                                                 #
#  NOMAD - Version 4 has been created and developed by                            #
#                 Viviane Rochon Montplaisir  - Polytechnique Montreal            #
#                 Christophe Tribes           - Polytechnique Montreal            #
#                                                                                 #
#  The copyright of NOMAD - version 4 is owned by                                 #
#                 Charles Audet               - Polytechnique Montreal            #
#                 Sebastien Le Digabel        - Polytechnique Montreal            #
#                 Viviane Rochon Montplaisir  - Polytechnique Montreal            #
#                 Christophe Tribes           - Polytechnique Montreal            #
#                                                                                 #
#  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             #
#  NSERC (Natural Sciences and Engineering Research Council of Canada),           #
#  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            #
#  for Data Valorization)                                                         #
#                                                                                 #
#  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     #
#  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       #
#  and Exxon Mobil.                                                               #
#                                                                                 #
#  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    #
#  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           #
#  Exxon Mobil.                                                                   #
#                                                                                 #
#  Contact information:                                                           #
#    Polytechnique Montreal - GERAD                                               #
#    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              #
#    e-mail: nomad@gerad.ca                                                       #
#                                                                                 #
#  This program is free software: you can redistribute it and/or modify it        #
#  under the terms of the GNU Lesser General Public License as published by       #
#  the Free Software Foundation, either version 3 of the License, or (at your     #
#  option) any later version.                                                     #
#                                                                                 #
#  This program is distributed in the hope that it will be useful, but WITHOUT    #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          #
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    #
#  for more details.                                                              #
#                                                                                 #
#  You should have received a copy of the GNU Lesser General Public License       #
#  along with this program. If not, see <http://www.gnu.org/licenses/>.           #
#                                                                                 #
#  You can find information on the NOMAD software at www.gerad.ca/nomad           #
#---------------------------------------------------------------------------------#

DESCRIPTION:

NOMAD is a C++ implementation of the Mesh Adaptive Direct Search (MADS)
algorithm, designed for constrained optimization of black-box functions.

The algorithms implemented are based on the book
"Derivative-Free and Blackbox Optimization", by Charles Audet and Warren Hare,
Springer 2017.


WEB PAGES:

https://www.gerad.ca/nomad/
Https://github.com/bbopt/nomad


CONTACT:

nomad@gerad.ca


VERSION WARNING:

This repository is for NOMAD 4. The previous version, NOMAD 3, 
is not on GitHub. NOMAD 3 will be deprecated in the future.

NOMAD 4 is similar in usage to NOMAD 3. 

NOMAD 4 has a new software architecture, uses OpenMP to run 
evaluations in parallel, and also has some new functionalities.

HOW TO CITE:

Please cite NOMAD 4 with reference:
C. Audet, S. Le Digabel, V. Rochon Montplaisir, and C. Tribes. 
Algorithm 1027: NOMAD version 4: Nonlinear optimization with the 
MADS algorithm. ACM Transactions on Mathematical Software
Volume 48, Issue 3, Article No.: 35, pp 1–22 
https://doi.org/10.1145/3544489


COMPILATION (Release):

On Linux, Unix, Windows and Mac OS X, NOMAD can be compiled using CMake.
The minimum version of CMake is 3.14. Older versions will trigger
an error. A recent C++ compiler is also required.

Before starting the procedure, we recommend setting the environment variable 
$NOMAD_HOME with the path where NOMAD has been copied. For Linux and OSX,

export NOMAD_HOME=/home/myUserName/PathToNomad

For Windows, add an environment variable %NOMAD_HOME% containing the path.
The remaining text uses the $NOMAD_HOME environment variable.

The procedure is the following. On the command line in the
 $NOMAD_HOME directory:

cmake -S . -B build/release     
    ---> Create the CMake files and directories for building (-B) in build/release.
         The source (-S) CMakeLists.txt file is in the $NOMAD_HOME directory.

         To deactivate compilation with OpenMP:
              cmake -DTEST_OPENMP=OFF -S . -B build/release

         To enable *C interface* building:
              cmake -DBUILD_INTERFACE_C=ON -S . -B build/release

         To enable *Matlab* interface building:
              cmake -DBUILD_INTERFACE_MATLAB=ON  -DTEST_OPENMP=OFF -S . -B build/release

              ! Before proceeding, have a look into 
              $NOMAD_HOME/interfaces/Matlab_MEX/readme.txt

              ! Extra flags might be required to prevent CMake errors. 

              ! Compiler version and Matlab version need to be compatible;
              ! Check https://www.mathworks.com/support/requirements/supported-compilers.html
               
              ! The Matlab interface will not be built if OpenMP is enabled.

              ! An extra addpath Matlab command must be done to have access 
              to nomad Mex binaries:
              addpath(strcat(getenv('NOMAD_HOME'),'/build/release/bin'))

         To enable *Python* interface (PyNomad) building:
              cmake -DBUILD_INTERFACE_PYTHON=ON -S . -B build/release

              ! Before proceeding, have a look into 
              $NOMAD_HOME/interfaces/PyNomad/readme.txt

              ! More details are provided in $NOMAD_HOME/interfaces/PyNomad/readme.txt 
             
              ! Building requires Cython. Cython can be obtained with
              Anaconda distribution platform.
  
              ! On *Windows*, using Visual Studio, see the user guide or the README
              to properly manage x86/x64 building of binaries. 
              
         To enable *Java* interface building (with Swig):
              cmake -DBUILD_INTERFACE_JAVA=ON -S . -B build/release

cmake --build build/release --config Release (for *Windows*)
or
cmake --build build/release (for *OSX* and *Linux*)
    ---> Build all the libraries and applications
    
         Option --parallel xx can be added for faster build.

         The option --config Release should be used on *Windows* 
         multi-configuration build environment (VisualStudio) to build only
         Release configuration. The default configuration is Debug.
         The same option should be used for *OSX* when using a *Xcode* project. 

cmake --install build/release --config Release (for *Windows*)
or
cmake --install build/release (for *OSX* and *Linux*)
    ---> Copy binaries and headers in build/release/[bin, include, lib]
         and in the examples/tests directories.

By default, the executable "nomad" will be installed into the directory:
build/release/bin/  (build/debug/bin/ when in debug mode). A symbolic link
is added in the bin directory (not functional for windows).

It is possible to build only a single application in its working directory:
(with NOMAD_HOME environment variable properly set)

cd $NOMAD_HOME/examples/basic/library/example1
cmake --build $NOMAD_HOME/build/release --target example1_lib.exe
cmake --install $NOMAD_HOME/build/release


COMPILATION (Debug):

The procedure to build the debug version is the following.
On the command line in the $NOMAD_HOME directory:

cmake -S . -B build/debug -D CMAKE_BUILD_TYPE=Debug
    ---> On *Windows*, all 4 configurations are configured
         (Debug, RelWithDebugInfo, MinSizeRel, Release); flag 
         CMAKE_BUILD_TYPE is ignored.

cmake --build build/debug --config Debug (for *Windows*)
or
cmake --build build/debug (for *OSX* and *Linux*)   
    ---> Build the libraries and applications
         
         Option --parallel xx can be added for faster build.

cmake --install build/debug --config Debug (for *Windows*)
or
cmake --install build/debug  (for *OSX* and *Linux*)
    ---> Copy binaries and headers in build/debug/[bin, include, lib]
         and in the examples/tests directories


EXAMPLES OF OPTIMIZATION:

Batch Mode:
There are examples in batch mode in $NOMAD_HOME/examples/basic/batch/.
In each directory, the blackbox functions (usually named bb) are compiled 
by default. The problem may be resolved using NOMAD and the parameter file:

$NOMAD_HOME/build/release/bin/nomad param.txt

For convenience, the path to $NOMAD_HOME/build/release/bin directory 
can be added to the $PATH environment variable. For *Windows*, this is 
achieved by setting the parameters for environment variable %PATH% to
%NOMAD_HOME\build\release\bin\
 

Library Mode:
There are examples in library mode in $NOMAD_HOME/examples/basic/library/.
In each directory, the executable may be compiled when building
the NOMAD application. The problems may be resolved by execution,
for instance:

./example_lib.exe

IMPORTANT: Library mode examples with *Windows* require to set 
the %PATH% environment variable (see above). Otherwise, the executables
cannot find NOMAD dlls.
