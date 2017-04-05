#!/usr/bin/env python
""" Setup script for the gwide package."""

__author__ = 'Tomasz W. Turowski'

from setuptools import setup, find_packages


setup(
    name='gwide',
    version='0.5.3',
    # py_modules=['gwide'],
    packages=find_packages(),
    install_requires=[
        'pypeaks',
        'pandas',
        'ruffus',
        'PyYAML',
        # 'matplotlib',
        # 'numpy'
        # 'plotly'
        # 'scipy'
    ],
    entry_points='''
        [console_scripts]
        gwidetRNA=gwide.gwidetRNA:tRNA
        gwidemRNA=gwide.gwidemRNA:mRNA
        gwiderRNA=gwide.gwiderRNA:rRNA
        gwideHittable=gwide.gwideHittable:hittable
        gwidePlot=gwide.gwidePlot:plot
        getFastaSeqs=gwide.parserTools:getFastaSeqs
        getGeneLength=gwide.parserTools:getGeneLength
        getIdFromName=gwide.parserTools:getIdFromName
        getNameFromId=gwide.parserTools:getNameFromId
        getGeneNamesFromGTF=gwide.parserTools:getGeneNamesFromGTF
        getNameFromId4Tab=gwide.parserTools:getNameFromId4Tab
    ''',
    scripts=[
        'gwide/scripts/novo2concat.py',
        'gwide/scripts/defineTerminator.py',
        'gwide/scripts/geneUsage.py',
        'gwide/scripts/codonCounter.py',
        'gwide/scripts/aminoacidCounter.py'
    ],
    author="Tomasz W. Turowski",
    description='Set of tools to downstream analysis of pyCRAC data',
    long_description='',
    author_email="twturowski@gmail.com",
    classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
)
