#!/usr/bin/env python
__author__ = 'tturowski'

#crating concat file from list of pileup files, adding experiment and project names as additional columns'''

import sys
import select
import time
import math
from optparse import OptionParser

def extract_names(line):
    path = line.strip()
    path_elem = path.split('/')
    if len(path_elem) > 1:
        exp_name = path_elem[len(path_elem)-2]
    else:
        exp_name = '.'
    return [exp_name, path]

def main():
    usage = "Usage: create pileups with pyPileup (pyCRAC package) then run one of below:"+"\n"+ \
            "A. When you have folder and subfolders, then subfolders names are taken as name of experiments \n" \
            "find . -name *.pileup | pileupsToConcat.py \n" \
            "or B. when is one folder: \n" \
            "find *.pileup | pileupsToConcat.py -e experiment_name" \
            ""
    parser = OptionParser(usage=usage)
    parser.add_option("-e", "--experiment_name", dest="exp_name", type="str",
                     help="experiment name", default=None)
    (options, args) = parser.parse_args()

    if not select.select([sys.stdin,],[],[],0.0)[0]: #expression from stackoverflow to check if there are data in standard input
        print usage
        print 'For more details use -h option.'
        exit()

    line_number = 0
    line_number_list = list()
    output_dict = dict()

    for line in sys.stdin:
        if not options.exp_name:
            names = extract_names(line)
            exp_name = str(names[0])
        else:
            try:
                exp_name = options.exp_name
            except:
                exp_name = '.'

        pileup_file = open(names[1])
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
    sys.stdout.write("# concat file from pileup files created: "+time.ctime()+"\n")
    sys.stdout.write("# gene\tposition\tnucleotide\thits\tsubstitutions\tdeletions\texperiment\thits_pM\tsubst_pM\tdel_pM\n")
    for i in line_number_list:
        sys.stdout.write('\t'.join(output_dict[str(i)]) + '\n')

main()
