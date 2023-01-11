#!/usr/bin/env python3
# 
# 2022 - 2023 Jan Provaznik (jan@provaznik.pro)
#
# Experimental support for statically linked PyNomad interface.

# Environment variables are used for configuration.
#
# (1) NOMAD_INSTALL for NOMAD installation directory.
# (2) NOMAD_VERSION for NOMAD version.
#
# Export the variables. 

import setuptools 
from Cython.Build import cythonize

import sys
import os.path
import os

nomad_install = os.environ['NOMAD_INSTALL']
nomad_version = os.environ['NOMAD_VERSION']

if not(nomad_install):
    print('Missing NOMAD_INSTALL env.')
    sys.exit(1)

if not(nomad_version):
    print('Missing NOMAD_VERSION env')
    sys.exit(1)

# Current support is limited to Linux (with GCC) only.

compile_args = []
compile_args.append('-w')
compile_args.append('-std=c++14')
compile_args.append('-pthread')

nomad_include_path = os.path.join(nomad_install, 'include/nomad')
nomad_library_path = os.path.join(nomad_install, 'lib')

link_args = []
link_args.append('-L{}'.format(nomad_library_path))
link_args.append('-l:libnomadAlgos.a')
link_args.append('-l:libnomadEval.a')
link_args.append('-l:libnomadUtils.a')
link_args.append('-l:libsgtelib.a')

setuptools.setup(
    name = 'PyNomad',
    version = nomad_version,
    author = 'Christophe Tribes',
    author_email = 'christophe.tribes@polymtl.ca',
    license = 'LGPL',
    description = 'Python interface to Nomad for blackbox optimization',
    url = 'https://github.com/bbot/nomad',
    ext_modules = cythonize(setuptools.Extension(
        'PyNomad',
        sources  =  [ 'PyNomad.pyx' ],
        libraries = [ 'gomp' ],
        include_dirs  =  [ nomad_include_path ],
        extra_compile_args  =  compile_args,
        extra_link_args  =  link_args,
        language  =  'c++'
    )),
    install_requires = [
        'Cython'
    ]
)

