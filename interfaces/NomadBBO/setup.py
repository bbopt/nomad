#!/usr/bin/env python3
#
# Support for statically linked NomadBBO interface.

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

import numpy


# Environment processing

env_nomad_src = os.environ.get('NOMAD_SRC')
env_nomad_build_dir = os.environ.get('NOMAD_BUILD_DIR')
env_nomad_msvc_flag = os.environ.get('NOMAD_MSVC_FLAG')
env_nomad_msvc_conf = os.environ.get('NOMAD_MSVC_CONF')
env_nomad_openmp_flag = (os.environ.get('BUILD_OPENMP') == 'TRUE')
env_nomad_sysroot = os.environ.get('SYSROOT')
env_nomad_archflags = os.environ.get('ARCHFLAGS')


if env_nomad_archflags:
    print('Building NomadBBO with ARCHFLAGS: ',os.environ.get('ARCHFLAGS'))

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

if env_nomad_openmp_flag:
    print('Building NomadBBO with OpenMP support.')

if env_nomad_sysroot:
    print('Building NomadBBO with isysroot set.')

# Construct base paths

# Create the list of include files for Nomad
path_include = [ env_nomad_src ]
path_include.append(numpy.get_include())


path_library_nomad = os.path.join(env_nomad_build_dir, 'src')
path_library_sgtelib = os.path.join(env_nomad_build_dir, 'ext', 'sgtelib')

# Compiler and linker configuration

setup_compile_args = []
setup_compile_args.append('-DNOMAD_STATIC_BUILD')
setup_link_args = []

if env_nomad_sysroot:
    clean_sysroot = env_nomad_sysroot.strip('"').strip("'")
    setup_compile_args.append('-isysroot')
    setup_compile_args.append(clean_sysroot)
    setup_link_args.append('-isysroot')
    setup_link_args.append(clean_sysroot)

if env_nomad_openmp_flag:
    setup_compile_args.append('-fopenmp')
    setup_link_args.append('-lgomp')

if env_nomad_msvc_flag:
    setup_compile_args.append('/std:c++17')
else:
    setup_compile_args.append('-w')
    setup_compile_args.append('-std=c++17')
    setup_compile_args.append('-pthread') 


# MSVC linker automagically resolves static libraries
# by their base names if given appropriate search paths.
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


path_version = os.path.join(env_nomad_src, 'nomad_version.hpp')

with open(path_version, encoding = 'utf-8') as file:
    pattern = '#define\s+NOMAD_VERSION_NUMBER\s+"([^"]+)"'
    for line in file:
        if (match := re.match(pattern, line)):
            nomad_version = match.group(1)
            break

# Off we go.

setuptools.setup(
    name="NomadBBO",
    author="Christophe Tribes and Edward Hallé-Hannan",
    author_email="christophe.tribes@polymtl.ca",
    description="Python interface to NOMAD (version " + nomad_version + ") for blackbox optimization with categorical variables",
    url="https://github.com/bbopt/nomad",

    packages=[
        "NomadBBO",
        "NomadBBO.python_src",
        "NomadBBO.python_src.data_managers",
        "NomadBBO.python_src.optimizers",
    ],
    package_dir={
        "NomadBBO": ".",
        "NomadBBO.python_src": "python_src",
        "NomadBBO.python_src.data_managers": "python_src/data_managers",
        "NomadBBO.python_src.optimizers": "python_src/optimizers",
    },
    include_package_data=True,

    ext_modules=cythonize([
        setuptools.Extension(
            "NomadBBO.NomadBBO_src",
            sources=["NomadBBO_src.pyx"],
            include_dirs=path_include,
            extra_compile_args=setup_compile_args,
            extra_link_args=setup_link_args,
            extra_objects=setup_extra_objects,
            libraries=setup_libraries,
            library_dirs=setup_library_dirs,
            language="c++",
        ),
    ]),
)
