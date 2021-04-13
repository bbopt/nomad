from distutils.core import setup, Extension
from Cython.Build import cythonize

import os
import sys


# Todo : investigate why it is not running properly when using OpenMP.

# Windows not supported
if sys.platform.startswith("win"):
    print("The", sys.platform, "platform is not supported.")
    exit()
    
if (len(sys.argv) != 5 and len(sys.argv) != 3):
    print("The script ", str(sys.argv[0]), " requires 4 arguments (building in place) or 2 arguments (installing). When building in place, arguments 1 and 2 are for passing Nomad options.")
    exit()

use_openmp=0
root_build_dir=""
os_include_dirs=""
if (len(sys.argv) == 5):
    #print("original sys argv: " + str(sys.argv))
    use_openmp = (int)(sys.argv[1])     # Argument 1 is the flag for using openMP or not
    root_build_dir = str(sys.argv[2])   # Argument 2 is to path to find libraries and headers
    os_include_dirs = [root_build_dir + "/../../src"]
    del sys.argv[1]
    del sys.argv[1]
    # print("new sys argv: " + str(sys.argv))

build_lib_dir = root_build_dir + "/src"
installed_lib_dir = root_build_dir + "/lib"

compile_args = []
compile_args.append("-std=c++14")
compile_args.append("-Wall")
compile_args.append("-Wextra")
compile_args.append("-pthread")

link_args = []
if (use_openmp==1):
    compile_args.append("-fopenmp")
    compile_args.append("-DUSE_OMP")
    link_args.append("-fopenmp")
    # print(compile_args)

# Use gcc
if (str(os.environ.get("CC")) == ""):
   os.environ["CC"] = "gcc"
if (str(os.environ.get("CXX")) == ""):
   os.environ["CXX"] = "g++"

# Look for librairies in Nomad distribution
# Default variables are for Linux
lib_extension="so"

# The rpath is set to the libraries installation directory.
link_args.append("-Wl,-rpath," + installed_lib_dir)

if sys.platform == "darwin":
     link_args.append("-headerpad_max_install_names")
     lib_extension = "dylib"

# For building/linking the interface, the libraries build directory is used.
link_args.append(build_lib_dir + "/libnomadUtils." + lib_extension)
link_args.append(build_lib_dir + "/libnomadEval." + lib_extension)
link_args.append(build_lib_dir+ "/libnomadAlgos." + lib_extension)

setup(
    name='PyNomad',
    ext_modules = cythonize(Extension(
           "PyNomad", # extension name
           sources = ["PyNomad.pyx"], # Cython source and interface
           include_dirs = os_include_dirs,
           extra_compile_args = compile_args,
           extra_link_args = link_args,
           language = "c++"
           ))
)

