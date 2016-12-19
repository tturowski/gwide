#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "1.0"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

from Bio import  SeqIO
from optparse import OptionParser
import argparse
from argparse import RawTextHelpFormatter
import select, os, re
import pandas as pd

usage = "Usage: calculate codon usage of fasta sequences"

parser = argparse.ArgumentParser(usage=usage, formatter_class=RawTextHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="fasta_file", help="Provide the path to your fasta file",
                 metavar="FILE", default=None)
args = parser.parse_args()

in_seq_handle = open(args)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()
seq_dict_keys =  seq_dict.keys()

print seq_dict_keys
exit()

a = 'ATGTATTATTAA'

codons = (a[n:n+3] for n in xrange(0,len(a),3)) # creates generator

dict_codons = {}

for codon in codons:
    if dict_codons.has_key(codon):
        dict_codons[codon] += 1
    else:
        dict_codons[codon] = 1

print dict_codons