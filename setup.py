from setuptools import setup

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: Implementation :: PyPy",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name='wdl_parser',
    version='0.0.1',
    description='Utilities for parsing .wdl files',
    url='https://github.com/ambrosejcarr/parse_wdl.git',
    author='Ambrose J. Carr',
    author_email='mail@ambrosejcarr.com',
    package_dir={'': 'src'},
    packages=['wdl_parser/test'],
    install_requires=['pyparsing'],
    classifiers=CLASSIFIERS,
    # include_package_data=True
)
