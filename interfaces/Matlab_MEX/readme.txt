The Matlab MEX interface allows to run NOMAD within the command line of Matlab.

Creating the Matlab MEX Interface to NOMAD requires to build source codes.
Building the interface requires compatibility of the versions of Matlab and
 the compiler that you want to use. Check the compatibility at:

https://www.mathworks.com/support/requirements/supported-compilers.html

************************
KNOWN ISSUES
************************

On Apple OSX computer with ARM64 and X86 architectures.
There might be incompatibility between the Matlab API binaries and
the binaries obtained when building Nomad. When build Nomad Mex interface 
this can result in a message like :
Undefined symbols for architecture arm64:
This issue can be resolved by forcing the architecture when configuring 
the build with a flag like -DCMAKE_OSX_ARCHITECTURES=x86_64.

************************
HOW TO BUILD AND INSTALL
************************
The interface build is managed by CMake that can be run *at NOMAD root*. 

The configuration command:
   cmake -DTEST_OPENMP=OFF -DBUILD_INTERFACE_MATLAB=ON -S . -B build/release.

Building the Matlab MEX interface is disabled when NOMAD uses OpenMP. 
Hence, the option -DTEST_OPENMP=OFF must be passed during configuration.

It may be required (Windows) to force the use of the 64bit version of the 
compiler with the command:
  cmake -DTEST_OPENMP=OFF -DBUILD_INTERFACE_MATLAB=ON -S . -B build/release -A x64

The command 
   cmake --build build/release 
or cmake --build build/release --config Release (for Windows) 
is used for building the selected configuration. 

The command
   cmake --install build/release 
must be run before using the Matlab nomadOpt function. 

Also, the Matlab command 
   addpath(strcat(getenv('NOMAD_HOME'),'/build/release/lib')) 
or addpath(strcat(getenv('NOMAD_HOME'),'/build/release/lib64')) 
must be executed to have access to the libraries and run the examples.

**********
HOW TO USE
**********
Some tests are proposed in the directory to check that everything 
is up and running.

All functionalities of NOMAD are available by using the nomadOpt 
function in Matlab command line. 

NOMAD parameters are provided in a Matlab structure with keywords and 
values using the same syntax as used in NOMAD parameter files. 
For example, 
   params = struct('initial_mesh_size','* 10','MAX_BB_EVAL','100');

*******************
COMPATIBILITY ISSUE
*******************
Even if the building process works smoothly, when running the nomadOpt() command in Matlab you may obtain an error such as "Invalid MEX_file .... libstdc++.so.6: version `GLIBCXX_3.4.26' not found" (example obtained on linux-Ubuntu). This is an indication that the versions of Matlab and the compiler may not be compatible (See in the preamble the web site to check compatibility).



