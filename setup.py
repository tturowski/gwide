__author__ = 'tturowski'
#!/usr/bin/env python
""" Setup script for the gwide_toolkit package."""

from setuptools import setup

setup(
    name='gwide_toolkit',
    version='0.1',
    py_modules=['gwideToolkit'],
    install_requires=[
        'pandas',
        'ruffus',
        'PyYAML'
    ],
    entry_points='''
        [console_scripts]
        gwideToolkit=gwideToolkit:main
    ''',
    author="Tomasz Turowski",
    description='Set of downstream tools for pyCRAC',
    long_description='Automated analysis pipeline for CRAC data, along with additional utilities',
    author_email="twturowski@gmail.com",
    classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache License 2.0',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
)