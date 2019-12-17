###################################################################################
#                                                                                 #
#                                    README                                       #
#                                                                                 #
#---------------------------------------------------------------------------------#
#  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                #
#                                                                                 #
#  NOMAD - Version 4.0.0 has been created by                                      #
#                 Viviane Rochon Montplaisir  - Polytechnique Montreal            #
#                 Christophe Tribes           - Polytechnique Montreal            #
#                                                                                 #
#  The copyright of NOMAD - version 4.0.0 is owned by                             #
#                 Sebastien Le Digabel        - Polytechnique Montreal            #
#                 Viviane Rochon Montplaisir  - Polytechnique Montreal            #
#                 Christophe Tribes           - Polytechnique Montreal            #
#                                                                                 #
#  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    #
#  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     #
#  Electrique and IVADO (The Institute for Data Valorization)                     #
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
#    phone : 1-514-340-6053 #6928                                                 #
#    fax   : 1-514-340-5665                                                       #
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

DISCLAIMER :

This is a Beta version of NOMAD 4.0.0. It represents work in progress and
subject to change without notice.


DESCRIPTION :

NOMAD is a C++ implementation of the Mesh Adaptive Direct Search (MADS)
algorithm, designed for constrained optimization of black-box functions.

The algorithms implemented are based on the book
"Derivative-Free and Blackbox Optimization", by Charles Audet and Warren Hare,
Springer 2017.


WEB PAGE :

https://www.gerad.ca/nomad/


CONTACT :

nomad@gerad.ca


USAGE :

NOMAD 4 usage is similar to NOMAD 3, including batch and library modes. 


COMPILATION :

On Linux, Unix, and Mac OS X, NOMAD can be compiled using the makefile
located in the "src" directory.

Make options for Nomad 4:

If OpenMP is not available:
make NOOMP=1

To compile without using sgtelib:
make USE_SGTELIB=0

To compile in Debug mode:
make VARIANT=debug

The executable "nomad" will be produced to directory:
build/release/bin/  (build/debug/bin/ when in debug mode).


EXAMPLES :

Batch Mode:
There are examples in batch mode in examples/basic/batch/.
In each directory, the blackbox (usually named bb.exe) may be compiled using the
provided makefile.
The problem may be resolved using NOMAD and the parameter file param.txt:
nomad param.txt

Library Mode:
There are examples in library mode in examples/basic/library/.
In each directory, the executable may be compiled using the provided makefile
(using the same flags as for the NOMAD compilation). The problems may be resolved
by execution, for instance:
example_lib.exe
