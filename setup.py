#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("spindle_dynamics", ["spindle_dynamics.pyx"])]

setup(name='kt_simul',
      version='0.6',
      description='Kinetochore dynamics simulation',
      author='Guillaume Gay',
      author_email='gllm.gay@gmail.com',
      url='https://github.com/glyg/Kinetochore-segregation',
      cmdclass = {'build_ext': build_ext},
      packages = ['kt_simul'],
      package_data = {'kt_simul': ['data/*.xml',
                                    '*.pyx', '*.pyd',
                                    'data/images/*.svg',
                                    'data/images/*.png']},
      )
