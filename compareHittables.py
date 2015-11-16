#!/usr/bin/env python
import numpy as np
import sys, select, os, re, math, argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import pandas.tools.rplot as rplot
import gwideToolkit.methods as gtk

# usage: create hittables using pyReadCounter then
# i.e. find *hittable* | countHittable.py

usage = "Calculates correlation cooeficient between two or more datasets. Usage: create hittables using pyReadCounter then run script in the folder containing hittables"
parser = argparse.ArgumentParser(usage=usage, formatter_class=RawTextHelpFormatter)
files = parser.add_argument_group('Options')
files.add_argument("-g", dest="gtf_file", help="Provide the path to your gtf file.",
                 type=str, default=None)
universal = parser.add_argument_group('Universal options')
universal.add_argument("-n", dest="normalized", action="store_true", help="Use when you want to work on data normalized 'reads per Milion'. Not working if you reading data.tab Default: False", default=False)
universal.add_argument("-w", dest="whole_name", action="store_true", help="As defauls scripts takes 'a_b_c' from a_b_c_hittable_reads.txt as experiment name. Use this option if your file names do not suit to this pattern. Default: False", default=False)
universal.add_argument("-c", dest="gene_class", action="store_true", help="Calculate Pearson coefficient for different classes separately. Default: False", default=False)
universal.add_argument("-o", dest="output", choices=["p", "s", "k", "a"], help="Select from following options: p - Pearson (standard correlation coefficient); s - Spearman rank correlation; k - Kendall Tau correlation coefficient; a - all at once", default="p")
options = parser.parse_args()

gtf = gtk.getGTF(options.gtf_file)

print "# compareHittables.py is running..."

paths = list()
experiments = list()
no_of_reads = dict()
gene_type = str()
gene_type_reads = int()
genes_name = list()
normalizator = float()
classes = list()

def define_experiments(paths_in):
    for path in paths_in:
        paths.append(path.strip())
        file_path = path.strip().split('/')
        file_name = file_path[len(file_path)-1]
        if options.whole_name == False:
            name = "_".join(file_name.split('_')[0:3]) #take fist three elements of file name as experiment name
        elif options.whole_name == True:
            name = file_name.strip('_hittable_reads.txt')
        experiments.append(name)
        experiments.sort()

# if not select.select([sys.stdin,],[],[],0.0)[0]: #expression from stackoverflow to check if there are data in standard input
files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('_hittable_reads.txt')]
define_experiments(files)
# else:
#     define_experiments(sys.stdin)

data = pd.DataFrame(columns=[['gene_name', 'gene_id', 'type']+experiments]) # initialize Pandas DataFrame
print "# Currently reading GTF file: "+ gtf
for line in open(gtf, 'r'):
# for line in open('/home/tturowski/seq_references/Saccharomyces_cerevisiae.EF4.74.shortChNames_with_PolIII_transcripts_extended_slop_intergenic_sort.gtf', 'r'):
    if not line.startswith('#'):
        line_elements = line.strip().split('\t')
        type = str(line_elements[1])
        if type not in classes:
            classes.append(type)
        try:
            gene_name = re.search("gene_name\s\"(.*?)\"", str(line_elements[8])).group(1)
        except:
            gene_name = re.search("gene_id\s\"(.*?)\"", str(line_elements[8])).group(1) # when there is no gene name
        gene_id = re.search("gene_id\s\"(.*?)\"", str(line_elements[8])).group(1)
        if gene_name not in genes_name:
            genes_name.append(gene_name)
            gene_data = pd.DataFrame([[gene_name, gene_id, type]+([0]*len(experiments))], columns=(['gene_name', 'gene_id', 'type']+experiments))
            data = data.append(gene_data, ignore_index=True)
data = data.set_index(['gene_name'])

#create file with no of reads
no_of_reads_file = open("no_of_reads.table",'w')
no_of_reads_file.write("# experiment"+'\t'+"mapped_reads"+'\t'+"total_reads"+'\n')

for path in paths:
    hittable = open(path.strip())
    file_path = path.strip().split('/')
    file_name = file_path[len(file_path)-1]
    if options.whole_name == False:
            name = "_".join(file_name.split('_')[0:3]) #take fist three elements of file name as experiment name
    elif options.whole_name == True:
            name = file_name.strip('_hittable_reads.txt')
    no_of_reads[name] = dict()
    print "# Currently reading: "+file_name+" ..."
    for line in hittable:
        if line.startswith('# total number of reads') and not line.startswith('# total number of reads without'):
            total_reads = int(filter(str.isdigit, line)) # no of reads
            no_of_reads[name]['total_reads'] = total_reads
        if line.startswith('# total number of single reads'):
            total_mapped_reads = int(filter(str.isdigit, line)) # no of mapped reads
            no_of_reads[name]['total_mapped_reads'] = total_mapped_reads
            normalizator = 1000000.0/total_mapped_reads
        if not line.startswith('#'):
            line_elements = line.strip().split('\t')
            if len(line_elements) == 4:
                gene_name, hits = line_elements[0], int(line_elements[1])
                if options.normalized == True:
                    data.loc[gene_name, name] = int(math.ceil(float(hits*normalizator)))
                else:
                    data.loc[gene_name, name] = hits
    # print no_of_reads
    no_of_reads_file.write(name+'\t'+str(no_of_reads[name]['total_mapped_reads'])+'\t'+str(no_of_reads[name]['total_reads'])+'\n')
# print "# Creating output.tab file..."
# data.to_csv('data.tab', sep='\t')
no_of_reads_file.close()


corr_dict = {"p" : "pearson" , "k" : "kendall" , "s" : "spearman"}
if options.output == 'a':
    print "# Calculating all correlations..."
    for i in corr_dict:
        print "# Calculating correlations("+corr_dict[i]+")..."
        matrix = data.corr(method=corr_dict[i],min_periods=1)
        matrix.to_csv("genome_wide_correlation_"+corr_dict[i]+".table", sep='\t')

#calculate Pearson for different types
        if options.gene_class == True:
            for this_type in classes:
                new_data = data[data.type == this_type]
                matrix = new_data.corr(method=corr_dict[i],min_periods=1)
                matrix.to_csv(this_type+"_correlation_"+corr_dict[i]+".table", sep='\t')
else:
    print "# Calculating correlations("+corr_dict[options.output]+")..."
    matrix = data.corr(method=corr_dict[options.output],min_periods=1)
    matrix.to_csv("genome_wide_correlation_"+corr_dict[options.output]+".table", sep='\t')

    #calculate Pearson for different types
    if options.gene_class == True:
        for this_type in classes:
            new_data = data[data.type == this_type]
            matrix = new_data.corr(method=corr_dict[options.output],min_periods=1)
            matrix.to_csv(this_type+"_correlation_"+corr_dict[options.output]+".table", sep='\t')
print "Done."