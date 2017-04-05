#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "1.0"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

### not finished ### abandoned

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
                                              "or ar_b_20 = aromatic aminoacids below 20% within window. Option for position 1: positive, negative, charged, polar, hydrophobic, aromatic"
                                              "Options for position 2: a - above or equal. Position 3 is percent within sliding window i.e 20 = 2/10 or 3/15"
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

#reading fasta file
in_seq_handle = open(args.fasta_file)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta")) #dictionary with
in_seq_handle.close()
seq_dict_keys = seq_dict.keys()

#function to check -a or config.list syntax
def filter_parser(input_str):
    #checks
    if len(input_str.split("_")) != 3:
        exit("Filter parameter incorrect: "+input_str+". Check -h for help.")
    if input_str.split("_")[0] not in ['positive', 'negative', 'charged', 'polar', 'hydrophobic', 'aromatic']:
        exit("Filter parameter incorrect: "+input_str+". Check -h for help.")
    if input_str.split("_")[1] not in ['a']:
        exit("Filter parameter incorrect: " + input_str + ". Check -h for help.")
    if int(input_str.split("_")[2]) not in range(0,100):
        exit("Filter parameter incorrect: " + input_str + ". Check -h for help.")
    aa_groups = { 'positive'    :   ['R', 'H', 'K'], #basic
                  'negative'    :   ['D', 'E'], #acidic
                  'charged'     :   ['R', 'H', 'K', 'D', 'E'],
                  'polar'       :   ['S', 'T', 'N', 'Q'],
                  'hydrophobic' :   ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
                  'aromatic'    :   ['F', 'W', 'Y']
                  }
    return {'aa'        : (input_str.split("_")[0], aa_groups[input_str.split("_")[0]]), #returns tuple
            'comp'      : input_str.split("_")[1], #string
            'percent'   : int(input_str.split("_")[2]) #integer
            }

#choose filter input
filter_list = list()
if args.aa_type: filter_list = [filter_parser(args.aa_type)] #from -a
elif not args.config_list: exit("Filter config.list missing. Add file using -c config.list or input single filtering requirement using -a option.")
else:
    for filter in open(args.config_list): filter_list.append(filter_parser(filter.strip())) #from config.list

output_matrix = pd.DataFrame()

for name in seq_dict_keys:
    gene_name = id_to_gene[gene_id]
    a = str(seq_dict[name].seq)
    windows = (a[n:n+args.window] for n in xrange(0,len(a)-args.window+1,1)) # creates generator
    no_of_windows = len(a)-args.window+1
    # print gene_name
    # print a

    for filter in filter_list:
        # print filter

        for window in enumerate(windows):
            # print window
            count = int()
            for aminoacid in filter['aa'][1]: count += window[1].count(aminoacid) #counting each aa from the given list in window
            # print count


    exit()
    #
    # if args.save_matrix == False:
    #     if args.id_given == True: #when gene ID are given as a names
    #         if dict_codons.has_key(args.codone): print id_to_gene[name]+'\t'+str(dict_codons[args.codone])+'\t'+args.codone
    #         else: print id_to_gene[name]+'\t'+str(0)+'\t'+args.codone
    #     elif args.id_given == False:
    #         if dict_codons.has_key(args.codone): print name+'\t'+str(dict_codons[args.codone])+'\t'+args.codone
    #         else: print name+'\t'+str(0)+'\t'+args.codone
    # else:
    #     dict_codons['sum']=len(dict_codons)
    #     matrix = pd.concat([matrix, pd.Series(data=dict_codons,name=name)], axis=1)
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