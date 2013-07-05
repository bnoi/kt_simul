# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

import kt_simul

install_requires = ['cython', 'numpy', 'numexpr',
                    'tables', 'pandas', 'matplotlib']

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

from Cython.Distutils import build_ext
from Cython.Distutils import Extension

ext_modules = [Extension("spindle_dynamics", ["spindle_dynamics.pyx"])]

setup(
    name='kt_simul',
    version=kt_simul.__version__,
    packages=find_packages(),
    author="BNOI Project",
    author_email="bnoi.project@gmail.com",
    description="""Python model of chromosome mouvements during mitosis in
                   Fission Yeast""",
    long_description=open('README.md').read(),
    install_requires=install_requires,
    include_package_data=True,
    url='https://github.com/bnoi/kt_simul.git',
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="CeCILL",
    entry_points={
        'console_scripts': [
            #'proclame-sm = sm_lib.core:proclamer',
        ],
    },
    cmdclass = {'build_ext': build_ext}
)
