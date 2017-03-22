#!/usr/bin/env python
__author__ = 'Hywel and Tomasz'
__copyright__	= "Copyright 2015"
__version__		= "0.2"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

from ruffus import *
import subprocess, os, argparse, time
import yaml

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
#seting up option parser
parser = argparse.ArgumentParser(description='Usage: ruffus scirpt to calculate Pearson correlation coeficient directly from *.novo files. Make new folder, cp or ln into all novofiles and run ruffus script. IMPORTANT: name of novo file should be name of experiments')
parser.add_argument("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     type=str, default=None)
parser.add_argument("-w", dest="whole_name", action="store_true", help="As defauls scripts takes 'a_b_c' from a_b_c.novo as experiment name. Use this option if your file names do not suit to this pattern. Default: False", default=False)
args = parser.parse_args()

gtf = getGTF(args.gtf_file)

print "# Using GTF file: " + gtf

#listing novo files
files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.novo')]  #gives list of files in current directory

#seting up no. of processes to use
number_of_processors_to_use = len(files) # or however many processors you want...

#running pyReadCounter on novofiles
@subdivide(files, suffix('.novo'), '*hittable_reads.txt')
def run_pyReadCounter(infile, outfile):
    print infile
    subprocess.call('pyReadCounters.py -f '+infile+' --gtf='+gtf, shell=True)
    subprocess.call(r'rm *count_output_reads.gtf *file_statistics_reads.txt *intron_and_UTR_overlap_reads.gtf', shell=True)

#creating tables
@merge(run_pyReadCounter, '*.table')
def merge_files(infiles, outfile):
    if args.whole_name == True:
        subprocess.call('compareHittables.py -c -w -g ' + gtf, shell=True)
    else:
        subprocess.call('compareHittables.py -c -g ' + gtf, shell=True)

pipeline_run(multiprocess=number_of_processors_to_use)

# printing summary
log = open("Pearson_log.txt", 'w')
log.write("Ruffus script run on "+time.strftime("%d/%m/%Y")+' at '+time.strftime("%H:%M:%S"))
log.write("Using GTF file: " + gtf)
log.close()