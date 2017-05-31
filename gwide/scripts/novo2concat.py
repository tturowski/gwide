#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "0.1"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

"""script using ruffus pipeline to do for each novofile in the folder: mkdir, and run pypileup,
then run pileupsToConcat.py and run script with plotting parameters: analyse_tRNA_from_concat.py"""

from ruffus import *
import subprocess, re, os, argparse, time, yaml, math, shutil
import gwide.methods as gtm

# def novo2concat():
#seting up option parser
parser = argparse.ArgumentParser(description='Usage: ruffus scirpt designed to make concat file from *.novo files. Make new folder, cp or ln into all novofiles and run novo2concat. IMPORTANT: name of novo file should be name of experiment')
parser.add_argument("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     type=str, default=None)
parser.add_argument("-t", "--tab_file", dest="tab_file", help="Provide the path to your tab genome file.",
                     type=str, default=None)
parser.add_argument("-r", dest="ranges", help="Set up ranges for pyPileup. Default = 250", default=250)
parser.add_argument("--3end", dest="three_end", help="Use pyPileup option --3end to only report counts for the 3' end of the reads. Default = False",
                    action="store_true", default=False)
parser.add_argument("--5end", dest="five_end",
                    help="Use pyPileup option --5end to only report counts for the 5' end of the reads. Default = False",
                    action="store_true", default=False)
parser.add_argument("-l", dest="list_file", help="Provide the FULL path to your gene_names.list file.", type=str, default=None, required=True)
parser.add_argument("--tree", dest="tree", help="If you want to leave tree of catalogs including pilups within. Default = None.",
                     action="store_true", default=False)
parser.add_argument("--anti", dest="anti",  help="Create additional concat file with antisense reads Default = None.",
                    action="store_true", default=False)
parser.add_argument("-p", dest="prefix", help="Prefix for concat file name", type=str, default="")
args = parser.parse_args()

gtf, tab, ranges = gtm.getGTF(args.gtf_file), gtm.getTAB(args.tab_file), str(args.ranges)
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

#setting up concat and log files nameq
concat_name = args.prefix+"_r"+ranges
if args.three_end == True:
    concat_name = concat_name + "_3end"
elif args.five_end == True:
    concat_name = concat_name + "_5end"

#seting up no. of processes to use
number_of_processors_to_use = len(files)

#making pileups in separate directories
@subdivide(
    links,
    formatter(),
    os.path.join("{path[0]}", "*.pileup")
)
def create_pileup_files(input_file, output_files):
    os.chdir(os.path.dirname(input_file))
    command_part = ''
    if args.three_end: command_part = command_part+' --3end'
    if args.five_end: command_part = command_part+' --5end'
    if args.anti: command_part = command_part + ' --anti_sense'

    subprocess.call(r'pyPileup.py --gtf='+gtf+' --tab='+tab+' -g '+args.list_file+' -f ' + input_file + ' -r '+ranges + command_part, shell=True)
    os.chdir(root_dir)

#creating concat file
@merge(create_pileup_files, os.path.join(root_dir, concat_name))
def merge_files(infiles, concat):
    print "# Creating concat file..."
    #below pileupsToConcac.py script is implemented

    def extract_names(line):
        path = line.strip()
        path_elem = path.split('/')
        exp_name = path_elem[len(path_elem)-2]
        return str(exp_name)

    line_number = 0
    line_number_list = list()
    output_dict = dict()

    print infiles

    for file_path in infiles:
        exp_name = extract_names(file_path)
        if "sense-reads" in file_path:
            pileup_file = open(file_path)
            for line in pileup_file:
                if line.startswith('# total number of mapped reads:'):
                    line_elem_1 = line.strip().split('\t')
                    total_mapped_reads = float(line_elem_1[1])
                    normalizator = 1000000.0/total_mapped_reads
                if not line.startswith('#'):
                    line_elements = line.strip().split('\t')
                    if len(line_elements) > 1:
                        line_number = line_number + 1
                        line_number_list.append(line_number)
                        line_list = list()
                        for i in line_elements:
                            line_list.append(str(i))
                        line_list.append(exp_name)
                        line_list.append(str(int(math.ceil(float(line_elements[3])*normalizator))))
                        line_list.append(str(int(math.ceil(float(line_elements[4])*normalizator))))
                        line_list.append(str(int(math.ceil(float(line_elements[5])*normalizator))))
                        output_dict[str(line_number)] = line_list

    line_number_list.sort()
    output = open(concat+".concat", 'w')
    output.write("# concat file from pileup files created: "+time.ctime()+"\n")
    output.write("# gene\tposition\tnucleotide\thits\tsubstitutions\tdeletions\texperiment\thits_pM\tsubst_pM\tdel_pM\n")
    for i in line_number_list:
        output.write('\t'.join(output_dict[str(i)]) + '\n')
    output.close()

    if args.anti == True:
        line_number = 0
        line_number_list = list()
        output_dict = dict()

        for file_path in infiles:
            exp_name = extract_names(file_path)
            if "anti_sense-reads" in file_path:
                pileup_file = open(file_path)
                for line in pileup_file:
                    if line.startswith('# total number of mapped reads:'):
                        line_elem_1 = line.strip().split('\t')
                        total_mapped_reads = float(line_elem_1[1])
                        normalizator = 1000000.0 / total_mapped_reads
                    if not line.startswith('#'):
                        line_elements = line.strip().split('\t')
                        if len(line_elements) > 1:
                            line_number = line_number + 1
                            line_number_list.append(line_number)
                            line_list = list()
                            for i in line_elements:
                                line_list.append(str(i))
                            line_list.append(exp_name)
                            line_list.append(str(int(math.ceil(float(line_elements[3]) * normalizator))))
                            line_list.append(str(int(math.ceil(float(line_elements[4]) * normalizator))))
                            line_list.append(str(int(math.ceil(float(line_elements[5]) * normalizator))))
                            output_dict[str(line_number)] = line_list

        line_number_list.sort()
        output = open(concat+"_antisense_reads.concat", 'w')
        output.write("# concat file from pileup files created: " + time.ctime() + "\n")
        output.write(
            "# gene\tposition\tnucleotide\thits\tsubstitutions\tdeletions\texperiment\thits_pM\tsubst_pM\tdel_pM\n")
        for i in line_number_list:
            output.write('\t'.join(output_dict[str(i)]) + '\n')
        output.close()

#version using dataframe but was extreamly slow!
    # columns = ['gene','position','nucleotide','hits','substitutions','deletions','exp_name','n_hits','n_substitutions','n_deletions']
    # data = pd.DataFrame(columns=columns)
    # for file_path in infiles:
    #     exp_name = extract_names(file_path) #extracts folder_name from file path; folder_name = name from name.novo
    #     pileup_file = open(file_path)
    #     print file_path
    #     for line in pileup_file:
    #         if line.startswith('# total number of mapped reads:'):
    #             line_elem_1 = line.strip().split('\t')
    #             total_mapped_reads = float(line_elem_1[1])
    #             normalizator = 1000000.0/total_mapped_reads
    #         if not line.startswith('#'):
    #             line_elements = line.strip().split('\t')
    #             if len(line_elements) > 1:
    #                 row_data = pd.DataFrame(columns=columns)
    #                 row_data['gene'], row_data['position'], row_data['nucleotide'], row_data['hits'], row_data['substitutions'], row_data['deletions'] = [line_elements[0]], [line_elements[1]], [line_elements[2]], [line_elements[3]], [line_elements[4]], [line_elements[5]]
    #                 row_data['exp_name'] = [exp_name]
    #                 row_data['n_hits'] = [str(int(math.ceil(float(line_elements[3])*normalizator)))]
    #                 row_data['n_substitutions'] = [str(int(math.ceil(float(line_elements[4])*normalizator)))]
    #                 row_data['n_deletions'] = [str(int(math.ceil(float(line_elements[5])*normalizator)))]
    #                 data.append(row_data)
    # print data
    # data.to_csv(concat, sep='\t')

pipeline_run(multiprocess=number_of_processors_to_use)

#deleting directories
if args.tree == False:
    for d in directories:
        shutil.rmtree(d)

# printing summary
log = open(concat_name+"_log.txt", 'w')
log.write("Ruffus script run on "+time.strftime("%d/%m/%Y")+' at '+time.strftime("%H:%M:%S")+'\n')
log.write("Using GTF file: " + gtf+'\n')
log.write("Using TAB genome file: " + tab+'\n')
log.write("Using gene_names.list file: " +args.list_file+'\n')
log.write("Ranges set up to: "+ranges+'\n')
log.close()

