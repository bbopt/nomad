###################################################################################
#                                                                                 #
#                                    README                                       #
#                                                                                 #
#---------------------------------------------------------------------------------#
#  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                #
#                                                                                 #
#  NOMAD - Version 4 has been created by                                          #
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


WEB PAGE:

https://www.gerad.ca/nomad/


CONTACT:

nomad@gerad.ca


VERSION WARNING:

This repository is for NOMAD 4. NOMAD 3 is not on GitHub.

NOMAD 4 is similar in usage to NOMAD 3. It does not have all
functionalities from NOMAD 3 yet. 

NOMAD 4 has a new software architecture, uses OpenMP to run 
evaluations in parallel, and also has some new functionalities.


COMPILATION (Release):

On Linux, Unix, Windows and Mac OS X, NOMAD can be compiled using CMake.
The minimum version of CMake is 3.14. Older versions will trigger
an error. A recent C++ compiler is also required.

The procedure is the following. On the command line in the
 $NOMAD_HOME directory:

cmake -S . -B build/release     
    ---> Create the CMake files and directories for building (-B) in build/release.
         The source (-S) CMakeLists.txt file is in the $NOMAD_HOME directory.

         To enable time stats build:
              cmake -DTIME_STATS=ON -S . -B build/release

         To enable C interface building:
              cmake -DBUILD_INTERFACE_C=ON -S . -B build/release

         To enable Matlab interface building:
              cmake -DBUILD_INTERFACE_MATLAB=ON -S . -B build/release
              ! Compiler version and Matlab version need to be compatible;
              ! Check https://www.mathworks.com/support/requirements/supported-compilers.html
               
              ! The Matlab interface will not be built if OpenMP is enabled.

              ! An extra addpath Matlab command must be done to have access 
              to nomad Mex binaries, 

         To enable Python interface (PyNomad) building:
              cmake -DBUILD_INTERFACE_PYTHON=ON -S . -B build/release
              ! The Matlab interface will not be built if OpenMP is enabled.

              ! Building requires to have Cython. Cython can be obtained with
              Anaconda distribution platform.
  
              ! On *Windows*, using Visual Studio, see the user guide to properly
              manage X86/X64 building of binaries. 

         To deactivate compilation with OpenMP:
              cmake -DTEST_OPENMP=OFF -S . -B build/release


cmake --build build/release     
    ---> Build all the libraries and applications
    
         Option --parallel xx can be added for faster build.

         Option --config Release should be used on *Windows* to build only
         Release configuration. The default configuration is Debug.


cmake --install build/release   
    ---> Copy binaries and headers in build/release/[bin, include, lib]
         and in the examples/tests directories.

         Option --config Release should be used on *Windows* to install 
         Release configuration. The default configuration is Debug.

By default, the executable "nomad" will installed into the directory:
build/release/bin/  (build/debug/bin/ when in debug mode). A symbolic link
is added in the bin directory.

It is possible to build only a single application in its working directory:
(with NOMAD_HOME environment variable properly set)

cd $NOMAD_HOME/examples/basic/library/example1
cmake --build $NOMAD_HOME/build/release --target example1_lib.exe
cmake --install $NOMAD_HOME/build/release


COMPILATION (Debug):

The procedure to build the debug version is the following.
On the command line in the $NOMAD_HOME directory:

cmake -S . -B build/debug -D CMAKE_BUILD_TYPE=Debug
    ---> On Windows, all 4 configurations are always build
         (Debug, RelWithDebugInfo, MinSizeRel, Release); flag 
         CMAKE_BUILD_TYPE is ignored.

cmake --build build/debug     
    ---> Build the libraries and applications
         
         Option --parallel xx can be added for faster build.

         On *Windows*, the default configuration is Debug.

make --install build/debug   
    ---> Copy binaries and headers in build/debug/[bin, include, lib]
         and in the examples/tests directories


EXAMPLES OF OPTIMIZATION:

Batch Mode:
There are examples in batch mode in examples/basic/batch/.
In each directory, the blackbox functions (usually named bb) are compiled 
by default. The problem may be resolved using NOMAD and the parameter file:

nomad param.txt

Library Mode:
There are examples in library mode in examples/basic/library/.
In each directory, the executable may be compiled when building
Nomad application. The problems may be resolved by execution,
for instance:

example_lib.exe
