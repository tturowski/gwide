#!/usr/bin/env python
""" Setup script for the gwide package."""

__author__ = 'Tomasz W. Turowski'

from setuptools import setup, find_packages


setup(
    name='gwide',
    version='0.3.4',
    # py_modules=['gwide'],
    packages=find_packages(),
    install_requires=[
        'pypeaks',
        'pandas',
        'ruffus',
        'PyYAML',
        'matplotlib',
        'numpy'
    ],
    entry_points='''
        [console_scripts]
        novo2concat.py=gwide.scripts.novo2concat:novo2concat
        defineTerminator.py=gwide.scripts.defineTerminator:main
        geneUsage.py=gwide.scripts.geneUsage:main
        gwidetRNA=gwide.gwidetRNA:tRNA
        gwideHittable=gwide.gwideHittable:hittable
        gwidePlot=gwide.gwidePlot:plot
        getFastaSeqs=gwide.parserTools:getFastaSeqs
        getGeneLength=gwide.parserTools:getGeneLength
        getIdFromName=gwide.parserTools:getIdFromName
        getNameFromId=gwide.parserTools:getNameFromId
        getGeneNamesFromGTF=gwide.parserTools:getGeneNamesFromGTF
    ''',
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
