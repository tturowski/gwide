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
import pandas as pd

usage = "Script calculates aminoacid composition within sliding window for all given proteins" \
        "Usage: calculate codon usage of fasta sequences"

parser = argparse.ArgumentParser(usage=usage, formatter_class=RawTextHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="fasta_file", help="Provide the path to your fasta file containing protein sequences. Required.",
                 metavar="FILE", default=None, required=True)
files.add_argument("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                   metavar="FILE", default=None)
files.add_argument("-w", dest="window", help="size of sliding window. Default = 10",
                   type=int, default=10)
files.add_argument("-a", dest="aa_type", help="Type of filter applied i.e. b_a_70 = basic aminoacides above or equal 70% within window "
                                              "or ar_b_20 = aromatic aminoacids below 20% within window. Option for position 1: a - acidic; b - basic; ar - aromatic"
                                              "Options for position 2: a - above or equal; b - below; e - equal. Position 3 is percent within sliding window i.e 20 = 2/10 or 3/15"
                                              "Not used when -c used",
                   type=str, default=None)
files.add_argument("-c", dest="config_list", help="Config.list. Default=False",
                   metavar="FILE", default=None)
files.add_argument("--id", dest="id_given", help="gene ID given instead of gene names", action="store_true", default=False)
args = parser.parse_args()

#reading GTF file to GTF parser and creating id_to_gene list
gtf = GTF2.Parse_GTF()
gtf.read_GTF(gtm.getGTF(args.gtf_file))

id_to_gene = dict()
for gene_name in gtf.genes:
    gene_id = gtf.genes[gene_name]['gene_id']
    id_to_gene[gene_id] = gene_name

gene_name = id_to_gene[gene_id]

#reading fasta file
in_seq_handle = open(args.fasta_file)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()
seq_dict_keys = seq_dict.keys()

#function to check -a or config.list syntax
def filter_parser(input_str):
    if len(input_str.split("_")) != 3: exit("Filter parameter incorrect: "+input_str+". Check -h for help.")
    if input_str.split("_")[0] not in ['a', 'b', 'ar']: exit("Filter parameter incorrect: "+input_str+". Check -h for help.")
    if input_str.split("_")[1] not in ['a', 'b', 'e']: exit(
        "Filter parameter incorrect: " + input_str + ". Check -h for help.")
    if int(input_str.split("_")[2]) not in range(0,100): exit(
        "Filter parameter incorrect: " + input_str + ". Check -h for help.")
    return input_str.split("_")[0], input_str.split("_")[1], input_str.split("_")[2]


output_matrix = pd.DataFrame()


for name in seq_dict_keys:
    a = pd.Series(str(seq_dict[name].seq))
    print a
    exit()
#     a = str(seq_dict[name].seq)
#     codons = (a[n:n+3] for n in xrange(0,len(a),3)) # creates generator
#     dict_codons = {}
#     for codon in codons:
#         if dict_codons.has_key(codon):
#             dict_codons[codon] += 1
#         else:
#             dict_codons[codon] = 1
#     if args.save_matrix == False:
#         if args.id_given == True: #when gene ID are given as a names
#             if dict_codons.has_key(args.codone): print id_to_gene[name]+'\t'+str(dict_codons[args.codone])+'\t'+args.codone
#             else: print id_to_gene[name]+'\t'+str(0)+'\t'+args.codone
#         elif args.id_given == False:
#             if dict_codons.has_key(args.codone): print name+'\t'+str(dict_codons[args.codone])+'\t'+args.codone
#             else: print name+'\t'+str(0)+'\t'+args.codone
#     else:
#         dict_codons['sum']=len(dict_codons)
#         matrix = pd.concat([matrix, pd.Series(data=dict_codons,name=name)], axis=1)
#
# #saving matrix with all codones
# if args.save_matrix == True:
#     # save numbers for each codone
#     matrix = matrix.transpose()
#     matrix = matrix.fillna(0)
#     matrix.to_csv("codone_composition.tab", sep='\t')
#
#     #calculate and save ratio for each codone
#     codones_sum = matrix.pop('sum')
#     ratio_matrix = matrix.div(codones_sum, axis=0)
#     ratio_matrix.to_csv("codone_composition_ratio.tab", sep='\t')