from setuptools import setup, find_packages

setup(
    name='kt_simul',
    version="2.0",
    packages=find_packages(),
    author="BNOI Project",
    author_email="bnoi.project@gmail.com",
    description="Python model of chromosome mouvements during mitosis",
    long_description=open('README.md').read(),
    include_package_data=True,
    url='https://github.com/bnoi/kt_simul.git',
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="CeCILL",
)
