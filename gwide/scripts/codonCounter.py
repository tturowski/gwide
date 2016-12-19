#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "1.0"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

from Bio import  SeqIO
from pyCRAC.Parsers import GTF2
import gwide.methods as gtm
import argparse
from argparse import RawTextHelpFormatter

usage = "Usage: calculate codon usage of fasta sequences"

parser = argparse.ArgumentParser(usage=usage, formatter_class=RawTextHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="fasta_file", help="Provide the path to your fasta file",
                 metavar="FILE", default=None)
files.add_argument("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                   metavar="FILE", default=None)
files.add_argument("-c", dest="codone", help="codone that want to count",
                   type=str, default='CGA')
args = parser.parse_args()

gtf = GTF2.Parse_GTF()
gtf.read_GTF(gtm.getGTF(args.gtf_file))

id_to_gene = dict()
for gene_name in gtf.genes:
    gene_id = gtf.genes[gene_name]['gene_id']
    id_to_gene[gene_id] = gene_name

gene_name = id_to_gene[gene_id]



in_seq_handle = open(args.fasta_file)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()
seq_dict_keys =  seq_dict.keys()

for name in seq_dict_keys:
    a = str(seq_dict[name].seq)
    codons = (a[n:n+3] for n in xrange(0,len(a),3)) # creates generator
    dict_codons = {}
    for codon in codons:
        if dict_codons.has_key(codon):
            dict_codons[codon] += 1
        else:
            dict_codons[codon] = 1
    if dict_codons.has_key(args.codone): print id_to_gene[name]+'\t'+str(dict_codons[args.codone])+'\t'+args.codone
    else: print id_to_gene[name]+'\t'+str(0)+'\t'+args.codoneecho