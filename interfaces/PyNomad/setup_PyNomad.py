from distutils.core import setup, Extension
from Cython.Build import cythonize

#import numpy as np
import os
import sys

# This value is be set when calling make with OPENMP=1
# By default, use OpenMP. This will clash if NOMAD was not compiled with OpenMP.
use_openmp = 1
if (len(sys.argv) > 3):
    use_openmp = (int)(sys.argv[1])
    del sys.argv[1]
    print("new sys argv: " + str(sys.argv))

# Windows not supported
if sys.platform.startswith("win"):
    print("The", sys.platform, "platform is not supported.")
    exit()
else:
    if (str(os.environ.get("NOMAD_HOME")) == "" or str(os.environ.get("NOMAD_HOME")) == "None"):
        print("A NOMAD_HOME environment variable is needed for building Nomad for Python (PyNomad)")
        exit()
    os_include_dirs = [str(os.environ.get("NOMAD_HOME")) + "/src"]

compile_args = []
compile_args.append("-std=c++14")
compile_args.append("-Wall")
compile_args.append("-Wextra")
compile_args.append("-pthread")
link_args = []
link_args.append(str(os.environ.get("NOMAD_HOME")) + "/build/release/lib/libnomadUtils.4.0.0.so")
link_args.append(str(os.environ.get("NOMAD_HOME")) + "/build/release/lib/libnomadEval.4.0.0.so")
link_args.append(str(os.environ.get("NOMAD_HOME")) + "/build/release/lib/libnomadAlgos.4.0.0.so")

if (use_openmp):
    compile_args.append("-fopenmp")
    compile_args.append("-DUSE_OMP")
    link_args.append("-fopenmp")

# Look for librairies in Nomad distribution
if sys.platform.startswith("linux"):
    link_args.append("-Wl,-rpath," + str(os.environ.get("NOMAD_HOME")) + "/build/release/lib")

# (OSX) Prevent error message when changing the location of libnomadUtils.so,
# libnomadEval.so and libnomadAlgos.so
if sys.platform == "darwin":
     link_args.append("-headerpad_max_install_names")
     os.environ["CC"] = "gcc"
     os.environ["CXX"] = "g++"

setup(
    ext_modules = cythonize(Extension(
           "PyNomad", # extension name
           sources = ["PyNomad.pyx", "nomadCySimpleInterface.cpp"], # Cython source and interface
           include_dirs = os_include_dirs,
           extra_compile_args = compile_args,
           extra_link_args = link_args,
           language = "c++"
           ))
)

if sys.platform == "darwin":
    for PyNomadLib in glob.glob("PyNomad*.so"):
        os.system("install_name_tool -change libnomad.so " + os.environ.get("NOMAD_HOME") + "/lib/libnomad.so " + PyNomadLib)

