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


# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

import peak_detection

install_requires = ['numpy', 'scipy', 'pandas', 'cython', 'scikit-image']

def install_requirements(install_requires):
    """
    Install third party libs in right order.
    """
    import subprocess
    import pip

    for package in install_requires:
        try:
            __import__(package)
        except:
            pip.main(['install', package])

install_requirements(install_requires)

setup(
    name='peak_detection',
    version=peak_detection.__version__,
    packages=find_packages(),
    author="BNOI Project",
    author_email="bnoi.project@gmail.com",
    description="""Python implementation of the Gaussian peak detection
                   described in Segré et al. Nature Methods (2008).
                   See https://github.com/bnoi/peak_detection for details.""",
    long_description=open('README.md').read(),
    install_requires=install_requires,
    include_package_data=True,
    url='https://github.com/bnoi/peak_detection',
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Image Recognition",
    ],
    entry_points={
        'console_scripts': [
            #'proclame-sm = sm_lib.core:proclamer',
        ],
    },
)
