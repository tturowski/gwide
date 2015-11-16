#!/usr/bin/env python
__author__ = 'Hywel and Tomasz'
__copyright__	= "Copyright 2015"
__version__		= "0.2"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

"""script using ruffus to do for each novofile in the folder: mkdir, and run pypileup,
then run pileupsToConcat.py and run script with plotting parameters: analyse_tRNA_from_concat.py"""

from ruffus import *
import subprocess, re, os, argparse, time

#seting up option parser
parser = argparse.ArgumentParser(description='Usage: ruffus scirpts are designed to make plots directly from *.novo files. Make new folder, cp or ln into all novofiles and run ruffus script. IMPORTANT: name of novo file should be name of experiment')
parser.add_argument("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     metavar="FILE", default=os.environ['GTF_PATH'])
parser.add_argument("-t", "--tab_file", dest="tab_file", help="Provide the path to your tab genome file.",
                     metavar="FILE", default=os.environ['TAB_PATH'])
parser.add_argument("-r", dest="ranges", help="Set up ranges for pyPileups. Default = 250", default=250)

args = parser.parse_args()

gtf, tab, ranges = args.gtf_file, args.tab_file, str(args.ranges)
print "Using GTF file: " + gtf
print "Using TAB genome file: " + tab

#listing novo files
files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.novo')]  #gives list of files in current directory
directories = [re.sub(r'.novo$', '', d) for d in files]
links = []
root_dir = os.getcwd()

#making directories
for f, d in zip(files, directories):
    os.mkdir(d)
    os.chdir(d)
    subprocess.call('ln -s ../' + f + ' ' + f, shell=True)
    links.append(os.path.abspath('./'+f))
    os.chdir(root_dir)

#setting up concat and log files name
concat_name = "tRNA_r"+ranges

#seting up no. of processes to use
number_of_processors_to_use = len(files) # or however many processors you want...

#making pileups in separate directories
@subdivide(
    links,
    formatter(),
    os.path.join("{path[0]}", "*.pileup")
)
def create_pileup_files(input_file, output_files):
    print input_file
    os.chdir(os.path.dirname(input_file))
    # subprocess.call(r'pyPileup.py -g ~/seq_references/yeast_genomic_tRNA_genes.list -f ' + input_file + ' -r 250', shell=True)
    subprocess.call(r'pyPileup.py --gtf='+gtf+' --tab='+tab+' -g /homes2/tturowski/seq_references/yeast_genomic_tRNA_genes.list -f ' + input_file + ' -r '+ranges, shell=True)
    subprocess.call(r'rm anti*', shell=True)
    os.chdir(root_dir)

#creating concat file
@merge(create_pileup_files, os.path.join(root_dir, concat_name+".concat"))
def merge_files(infiles, tRNA_concat):
    python_path_filepath = os.path.join(root_dir, 'python_paths.txt')
    with open(python_path_filepath, 'w') as python_path_file:
        for line in infiles:
            python_path_file.write(line+'\n')
    subprocess.call('cat ' + python_path_filepath + ' | pileupsToConcat.py > ' + tRNA_concat, shell=True)

#plotting
@transform(merge_files, suffix('.concat'), '*.png')
def plot_raw(infile, outfile):
    print infile
    subprocess.call('analyse_tRNA_from_concat.py -c '+infile+' -g '+gtf+' -p raw '+' --5flank='+ranges+' --3flank='+ranges, shell=True)

@transform(merge_files, suffix('.concat'), '*.png')
def plot_normalized(infile, outfile):
    subprocess.call('analyse_tRNA_from_concat.py -c '+infile+' -g '+gtf+' -o fig -n -p normalized'+' --5flank='+ranges+' --3flank='+ranges, shell=True)

pipeline_run(multiprocess=number_of_processors_to_use)

# printing summary
log = open(concat_name+"_log.txt", 'w')
log.write("Ruffus script run on "+time.strftime("%d/%m/%Y")+' at '+time.strftime("%H:%M:%S"))
log.write("Using GTF file: " + gtf)
log.write("Using TAB genome file: " + tab)
log.write("Ranges set up to: "+ranges)
log.close()