# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

from setuptools import setup, find_packages

import kt_simul

setup(
    name='kt_simul',
    version="2.0",
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
    }
)
