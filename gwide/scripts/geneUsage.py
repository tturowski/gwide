#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "1.0"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

from tRNAFromConcatv2 import *
from optparse import OptionParser
import argparse
from argparse import RawTextHelpFormatter
import select, os, re
import pandas as pd
def main():
    usage = "Usage: calculate gene usage for tRNA genes"

    parser = argparse.ArgumentParser(usage=usage, formatter_class=RawTextHelpFormatter)
    files = parser.add_argument_group('Options for input files')
    files.add_argument("-w", dest="w_file", help="Provide the path to your tab file with anticodones wages",
                     metavar="FILE", default=None)
    files.add_argument("-e", dest="e_file", help="Provide the path to your tab file tRNA expression",
                     metavar="FILE", default=None)
    options = parser.parse_args()

    #function from Internet by Tommy Tang
    def ReverseComplement(seq):
        seq_dict = {'A':'U','T':'A','G':'C','C':'G', 'a':'u','t':'a','g':'c','c':'g'}
        return "".join([seq_dict[base] for base in reversed(seq)])

    if not options.w_file or not options.e_file:
        print "Some file(s) is missing, use -w and -e"
        exit()

    #making revere complement (DNA to RNA)
    wages = pd.read_csv(options.w_file, sep='\t')
    wages['wage'] = wages['wage'].convert_objects(convert_numeric=True)
    for i, c in enumerate(wages['codone']):
        wages.loc[i, 'anticodone'] = ReverseComplement(c)

    #calculating % of codon usage
    new_wages = pd.DataFrame()
    for isotype in sorted(set(wages['isotype'].tolist())):
        temp_df = wages[wages.isotype == isotype]
        temp_df['wage_percent'] = temp_df['wage'].divide(temp_df.sum()['wage'], axis=0, level=None, fill_value=None)
        temp_df['wage_percent_correction'] = 0
        new_wages = new_wages.append(temp_df)
    wages = new_wages

    #calculating average hits and getting codone and anticodone from gene ID
    expression = pd.read_csv(options.e_file, sep='\t')
    expression['average'] = expression.mean(axis=1, numeric_only=True)
    for i, gene_id in enumerate(expression['gene_id']):
        expression.loc[i, 'anticodone'] = gene_id[3:6]
        expression.loc[i, 'isotype'] = gene_id[1]

    #calculating % of gene usage for each anticodone
    new_expression = pd.DataFrame()
    for anticodone in set(expression['anticodone'].tolist()):
        temp_df = expression[expression.anticodone == anticodone]
        temp_df['exp_percent_anticodone'] = temp_df['average'].divide(temp_df.sum()['average'], axis=0, level=None, fill_value=None)
        new_expression = new_expression.append(temp_df)
    expression = new_expression

    # #calculating % of gene usage for each isotype
    # new_expression2 = pd.DataFrame()
    # for isotype in set(expression['isotype'].tolist()):
    #     temp_df = expression[expression.isotype == isotype]
    #     temp_df['exp_percent_isotype'] = temp_df['average'].divide(temp_df.sum()['average'], axis=0, level=None, fill_value=None)
    #     new_expression2 = new_expression2.append(temp_df)
    # expression = new_expression2

    #getting all alternative (wobble driven) codones in decoding
    a = wages[['anticodone','codone', 'isotype', 'wage_percent']]
    b = expression[['anticodone','isotype']]
    b['gene'] = 'gene'  #mark gene encoded anticodones
    c = pd.ordered_merge(a,b)
    d = c[c.gene != 'gene']     #leave only wobble-driven codones
    d = d.sort('isotype')

    #re-asigning codones to new anticodones and adding wages for some genes
    wobble_driven_anticodones = list()
    wobble_dict = {'A':'A','T':'G','G':'U','C':'A', 'a':'a','t':'G','g':'u','c':'a'}
    for i, row in d.iterrows():
        if row['isotype'] != 'Z':
            anti1 = wobble_dict[row['codone'][2]]
            anti23 = ReverseComplement(row['codone'])[1:3]
            d.loc[i, 'recog_ant'] = anti1+anti23
            wages.loc[wages.anticodone == anti1+anti23,'wage_percent_correction'] += row['wage_percent']
    wages['wages_sum'] = wages['wage_percent']+wages['wage_percent_correction']
    # print d
    d.to_csv('wages', sep='\t')
    # wages = new_wages.drop('isotype', axis=1)
    # expression = new_expression.drop('isotype', axis=1)

    olo =  pd.ordered_merge(expression,new_wages)
    # print olo
    olo.to_csv('output', sep='\t')

    print "Done."