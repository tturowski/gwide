#!/usr/bin/env python
__author__ = 'Hywel and Tomasz'
__copyright__	= "Copyright 2015"
__version__		= "0.2"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

"""script using ruffus to do for each novofile in the folder: mkdir, and run pypileup,
then run pileupsToConcat.py and run script with plotting parameters: analyse_rRNA_from_concat.py"""

from ruffus import *
import subprocess, re, os, argparse

parser = argparse.ArgumentParser(description='Usage: ruffus scirpts are designed to make plots directly from *.novo files. Make new folder, cp or ln into all novofiles and run ruffus script. IMPORTANT: name of novo file should be name of experiment')
parser.add_argument("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     metavar="FILE", default=os.environ['GTF_PATH'])
parser.add_argument("-t", "--tab_file", dest="tab_file", help="Provide the path to your tab genome file.",
                     metavar="FILE", default=os.environ['TAB_PATH'])
args = parser.parse_args()

gtf, tab = args.gtf_file, args.tab_file
print "Using GTF file: " + gtf
print "Using TAB genome file: " + tab

files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.novo')]  #gives list of files in current directory
directories = [re.sub(r'.novo$', '', d) for d in files]
links = []
root_dir = os.getcwd()
for f, d in zip(files, directories):
    os.mkdir(d)
    os.chdir(d)
    subprocess.call('ln -s ../' + f + ' ' + f, shell=True)
    links.append(os.path.abspath('./'+f))
    os.chdir(root_dir)

number_of_processors_to_use = len(files) # or however many processors you want...

@subdivide(
    links,
    formatter(),
    os.path.join("{path[0]}", "*.pileup")
)
def create_pileup_files(input_file, output_files):

    print input_file
    os.chdir(os.path.dirname(input_file))
    subprocess.call(r'pyPileup.py --gtf='+gtf+' --tab='+tab+' -g /homes2/tturowski/seq_references/rRNA_gene.list -f ' + input_file + ' -r 1000', shell=True)
    subprocess.call(r'rm anti*', shell=True)
    os.chdir(root_dir)

@merge(create_pileup_files, os.path.join(root_dir, 'rRNA.concat'))
def merge_files(infiles, rRNA_concat):
    python_path_filepath = os.path.join(root_dir, 'python_paths.txt')
    with open(python_path_filepath, 'w') as python_path_file:
        for line in infiles:
            python_path_file.write(line+'\n')
    subprocess.call('cat ' + python_path_filepath + ' | pileupsToConcat.py > ' + rRNA_concat, shell=True)

@transform(merge_files, suffix('.concat'), '*.png')
def plot_raw(infile, outfile):
    print infile
    subprocess.call('analyse_rRNA_from_concat.py -c '+infile+' -g '+gtf+' -p raw', shell=True)

@transform(merge_files, suffix('.concat'), '*.png')
def plot_normalized(infile, outfile):
    subprocess.call('analyse_rRNA_from_concat.py -c '+infile+' -g '+gtf+' -n -p normalized', shell=True)

pipeline_run(multiprocess=number_of_processors_to_use)