__author__ = 'tturowski'

import os, yaml, sys
import numpy as np
import pandas as pd

def getRefFile(file_from_options, file_type):
    """
    Sorting out source of gtf, fasta or tab path, in order (1) from options parser, (2) from ~/bin/default.aml
     or (3) from environmental variable $xxx_PATH
    :param file_from_options: file path from options parser, can be an empty string
                file_type: 'GTF', 'FASTA', 'TAB'
    :return: path to the file
    """
    file_types = {'GTF' : 'GTF_PATH', 'FASTA' : 'FASTA_PATH', 'TAB' : 'TAB_PATH'}
    if file_from_options:
        return file_from_options
    elif 'default.aml' in os.listdir(os.getenv("HOME")+'/bin/'):
        default = yaml.load(open(os.getenv("HOME")+'/bin/default.aml'))
        # print "# Using "+file_type+" file from ~/bin/default.aml"
        return default[file_types[file_type]]
    else:
        if file_types[file_type] in os.environ:
            # print "# Using "+file_type+" file from $GTF_PATH variable"
            return os.environ[file_types[file_type]]
        else:
            exit('Provide '+file_type+' file path using -g or setup default.aml in your bin folder')

def getGTF(gtf_from_options):
    return getRefFile(gtf_from_options, 'GTF')

def getFASTA(fasta_from_options):
    return getRefFile(fasta_from_options, 'FASTA')

def getTAB(tab_from_options):
    return getRefFile(tab_from_options, 'TAB')

def list_paths_in_current_dir(suffix=str(), stdin=False):
    """
    :param: suffix  -   lists paths in current directory ending with an indicated suffix only
            stdin   -   read from standard input instead current directory
    :return: list of paths in current dir ending with suffix
    """
    #choosing between curr dir and std input
    if stdin == False:
        where = os.listdir('.')
    elif stdin == True:
        where = sys.stdin
    paths = [f for f in where if os.path.isfile(f) and f.endswith(suffix)]
    return paths

def define_experiments(paths_in, whole_name=False):
    """
    Parse file names and extract experiment name from them
    :param whole_name As defaults script takes first 'a_b_c' from a_b_c_hittable_reads.txt as experiment name.
    :return: list of experiment names and list of paths.
    """
    experiments = list()
    paths = list()
    for path in paths_in:
        paths.append(path.strip())
        file_path = path.strip().split('/')
        file_name = file_path[len(file_path)-1]
        if whole_name == False:
            name = "_".join(file_name.split('_')[0:3]) #take fist three elements of file name as experiment name
        elif whole_name == True:
            name = file_name.strip('_hittable_reads.txt')
        experiments.append(name)
        if len(experiments) != len(paths):
            exit("No. of experiments is not equal to no. of paths")
    return experiments, paths