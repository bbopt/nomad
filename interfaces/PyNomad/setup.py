#!/usr/bin/env python3
#
# Support for statically linked PyNomad interface.

# Environment variables are used for configuration.
#
# (1) NOMAD_SRC for NOMAD source directory.
# (2) NOMAD_BUILD for NOMAD build directory.
# (3) NOMAD_MSVC_FLAG for signalization that MSVC will be used.
# (4) NOMAD_MSVC_CONF for signalization of MSVC build configuration.
#
# Export the variables. 

import setuptools 
from Cython.Build import cythonize

import sys
import os.path
import os
import re

# Environment processing

env_nomad_src = os.environ.get('NOMAD_SRC')
env_nomad_build_dir = os.environ.get('NOMAD_BUILD_DIR')
env_nomad_msvc_flag = os.environ.get('NOMAD_MSVC_FLAG')
env_nomad_msvc_conf = os.environ.get('NOMAD_MSVC_CONF')

if not(env_nomad_src):
    print('Missing NOMAD_SRC env.')
    sys.exit(1)
    
if not(env_nomad_build_dir):
    print('Missing NOMAD_BUILD_DIR env.')
    sys.exit(1)

if not(env_nomad_msvc_flag):
    print('Missing NOMAD_MSVC_FLAG, assuming GCC compatible compiler.')

if not(env_nomad_msvc_conf):
    print('Missing NOMAD_MSVC_CONF, assuming Release configuration.')

# Construct base paths

path_include = env_nomad_src
path_library_nomad = os.path.join(env_nomad_build_dir, 'src') 
path_library_sgtelib = os.path.join(env_nomad_build_dir, 'ext', 'sgtelib') 

# Compiler and linker configuration

setup_compile_args = []
setup_compile_args.append('-DNOMAD_STATIC_BUILD')

if env_nomad_msvc_flag:
    setup_compile_args.append('/std:c++14')
else:
    setup_compile_args.append('-w')
    setup_compile_args.append('-std=c++14')
    setup_compile_args.append('-pthread')

    # Use GCC on non-MSVC builds.
    if (os.environ.get("CC") == ""):
       os.environ["CC"] = "gcc"
    if (os.environ.get("CXX") == ""):
       os.environ["CXX"] = "g++"

# MSVC linker automagically resolves static libraries
# by their base names if given apropriate search paths.
#
# Hence setup_libraries and setup_library_dirs.

setup_libraries = []
setup_library_dirs = []

# GCC compatible linkers can be explicitly instructed to use 
# static libraries as well (-l:libname.a). It is also possible
# to provide the static libraries directly as extra objects
# for the linker.
# 
# Hence setup_extra_objects.

setup_extra_objects = []

if env_nomad_msvc_flag:
    setup_libraries.append('nomadStatic')
    setup_libraries.append('sgtelibStatic')

    setup_library_dirs.append(os.path.join(path_library_nomad, env_nomad_msvc_conf))
    setup_library_dirs.append(os.path.join(path_library_sgtelib, env_nomad_msvc_conf))

else:
    setup_extra_objects.append(os.path.join(path_library_nomad, 'libnomadStatic.a'))
    setup_extra_objects.append(os.path.join(path_library_sgtelib, 'libsgtelibStatic.a'))

# Determine package version directly from NOMAD header file.

__version__ = "unknown"

path_version = os.path.join(path_include, 'nomad_version.hpp')

with open(path_version, encoding = 'utf-8') as file:
    pattern = '#define\s+NOMAD_VERSION_NUMBER\s+"([^"]+)"'
    for line in file:
        if (match := re.match(pattern, line)):
            __version__ = match.group(1)
            break

# Off we go.

setuptools.setup(
    name = 'PyNomad',
    version = __version__,
    author = 'Christophe Tribes',
    author_email = 'christophe.tribes@polymtl.ca',
    license = 'LGPL',
    description = 'Python interface to Nomad for blackbox optimization',
    url = 'https://github.com/bbopt/nomad',
    ext_modules = cythonize(setuptools.Extension(
        'PyNomad',
        sources = [ 'PyNomad.pyx' ],
        include_dirs = [ path_include ],
        extra_compile_args = setup_compile_args,
        extra_objects = setup_extra_objects,
        libraries = setup_libraries,
        library_dirs = setup_library_dirs,
        language  =  'c++'
    )),
)

