#!/usr/bin/env python3
# 
# 2022 - 2023 Jan Provaznik (jan@provaznik.pro)
#
# Experimental support for statically linked PyNomad interface.

# Environment variables are used for configuration.
#
# (1) NOMAD_INSTALL for NOMAD installation directory.
# (2) NOMAD_VERSION for NOMAD version.
# (3) NOMAD_VIAMSVC for signalization that MSVC will be used.
#
# Export the variables. 

import setuptools 
from Cython.Build import cythonize

import sys
import os.path
import os

# Environment processing

env_nomad_install = os.environ.get('NOMAD_INSTALL')
env_nomad_version = os.environ.get('NOMAD_VERSION')
env_nomad_viamsvc = os.environ.get('NOMAD_VIAMSVC')

if not(env_nomad_install):
    print('Missing NOMAD_INSTALL env.')
    sys.exit(1)

if not(env_nomad_version):
    print('Missing NOMAD_VERSION env')
    sys.exit(1)

if not(env_nomad_viamsvc):
    print('Missing NOMAD_VIAMSVC, assuming GCC compatible compiler.')

# Construct paths

nomad_include_path = os.path.join(env_nomad_install, 'include', 'nomad')
nomad_library_path = os.path.join(env_nomad_install, 'lib')

nomad_library_list = [ 
    'nomadStatic', 
    'sgtelibStatic' 
]

# Compiler and linker configuration

setup_compile_args = []
setup_compile_args.append('-DNOMAD_STATIC_BUILD')

if env_nomad_viamsvc:
    setup_compile_args.append('/std:c++14')
else:
    setup_compile_args.append('-w')
    setup_compile_args.append('-std=c++14')
    setup_compile_args.append('-pthread')

setup_libraries = []
setup_library_dirs = []
setup_extra_objects = []

if env_nomad_viamsvc:
    setup_libraries.extend(nomad_library_list)
    setup_library_dirs.append(nomad_library_path)
else:
    for library_name in nomad_library_list:
        library_path = os.path.join(nomad_library_path, f'lib{library_name}.a') 
        setup_extra_objects.append(library_path)

# Off we go.

setuptools.setup(
    name = 'PyNomad',
    version = env_nomad_version,
    author = 'Christophe Tribes',
    author_email = 'christophe.tribes@polymtl.ca',
    license = 'LGPL',
    description = 'Python interface to Nomad for blackbox optimization',
    url = 'https://github.com/bbot/nomad',
    ext_modules = cythonize(setuptools.Extension(
        'PyNomad',
        sources = [ 'PyNomad.pyx' ],
        include_dirs = [ nomad_include_path ],
        extra_compile_args = setup_compile_args,
        extra_objects = setup_extra_objects,
        libraries = setup_libraries,
        library_dirs = setup_library_dirs,
        language  =  'c++'
    )),
)

