#!/usr/bin/env python
import numpy as np
import sys, select, os, re, math, argparse, yaml
from argparse import RawTextHelpFormatter
import pandas as pd

# usage: create hittables using pyReadCounter then
# i.e. find *hittable* | countHittable.py

#sorting out GTF file
def getGTF(gtf_from_options):
    if gtf_from_options:
        return gtf_from_options
    elif 'default.aml' in os.listdir(os.getenv("HOME")+'/bin/'):
        default = yaml.load(open(os.getenv("HOME")+'/bin/default.aml'))
        # print "# Using GTF file from ~/bin/default.aml"
        return default['GTF_PATH']
    else:
        if os.environ['GTF_PATH']:
            # print "# Using GTF file from $GTF_PATH variable"
            return os.environ['GTF_PATH']
        else:
            exit('Provide GTF file path using -g or setup default.aml in your bin folder')

#setup option parser
usage = "Usage: prints genom wide plot, for genes 'as you wish' and aligned 'as you wish'"
parser = argparse.ArgumentParser(usage=usage, formatter_class=RawTextHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                 type=str, default=None)
options = parser.parse_args()

gtf_file = getGTF(options.gtf_file)

paths = list()
experiments = list()
no_of_reads = dict()
gene_type = str()
gene_type_reads = int()
genes_name = list()
normalizator = float()

for path in sys.stdin:
    paths.append(path)
    file_path = path.strip().split('/')
    file_name = file_path[len(file_path)-1]
    name = "_".join(file_name.split('_')[0:3]) #take fist three elements of file name as experiment name
    experiments.append(name)

data = pd.DataFrame(columns=[['gene_name', 'gene_id', 'type']+experiments]) # initialize Pandas DataFrame

print "Currently reading GTF file: "+os.environ['GTF_PATH']
for line in open(os.environ['GTF_PATH'], 'r'):
# for line in open('/home/tturowski/seq_references/Saccharomyces_cerevisiae.EF4.74.shortChNames_with_PolIII_transcripts_extended_slop_intergenic_sort.gtf', 'r'):
    if not line.startswith('#'):
        line_elements = line.strip().split('\t')
        type = str(line_elements[1])
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
data.to_csv('gtf_read_TEMP.tab', sep='\t')

for path in paths:
    hittable = open(path.strip())
    file_path = path.strip().split('/')
    file_name = file_path[len(file_path)-1]
    name = "_".join(file_name.split('_')[0:3]) #take fist three elements of file name as experiment name
    print "Currently reading: "+file_name
    for line in hittable:
        if line.startswith('# total number of single reads'):
            total_mapped_reads = int(filter(str.isdigit, line)) # no of mapped reads
            no_of_reads[name] = total_mapped_reads
            normalizator = 1000000.0/total_mapped_reads
        # if line.startswith('##'):
        #     gene_type = line.strip("#").strip().split('\t')[0]
        #     gene_type_reads = int(line.strip("#").strip().split('\t')[1])
        if not line.startswith('#'):
            line_elements = line.strip().split('\t')
            if len(line_elements) == 4:
                gene_name, hits = line_elements[0], int(line_elements[1])
                # print gene_name
                data.loc[gene_name, name] = int(math.ceil(float(hits*normalizator)))
print "Creating output.tab file..."
data.to_csv('output.tab', sep='\t')