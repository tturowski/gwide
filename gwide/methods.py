__author__ = 'tturowski'

import os, yaml, sys
import numpy as np
import pandas as pd

def getGTF(gtf_from_options):
    """
    Sorting out source of gtf path, in order (1) from options parser, (2) from ~/bin/default.aml
     or (3) from environmental variable $GTF_PATH
    :param gtf_from_options: gtf path from options parser, can be an empty string
    :return: path to gtf file
    """
    if gtf_from_options:
        return gtf_from_options
    elif 'default.aml' in os.listdir(os.getenv("HOME")+'/bin/'):
        default = yaml.load(open(os.getenv("HOME")+'/bin/default.aml'))
        # print "# Using GTF file from ~/bin/default.aml"
        return default['GTF_PATH']
    else:
        if 'GTF_PATH' in os.environ:
            # print "# Using GTF file from $GTF_PATH variable"
            return os.environ['GTF_PATH']
        else:
            exit('Provide GTF file path using -g or setup default.aml in your bin folder')

def getFASTA(fasta_from_options):
    """
    Sorting out source of fasta path, in order (1) from options parser, (2) from ~/bin/default.aml
     or (3) from environmental variable $FASTA_PATH
    :param fasta_from_options: fasta path from options parser, can be an empty string
    :return: path to fasta file
    """
    if fasta_from_options:
        return fasta_from_options
    elif 'default.aml' in os.listdir(os.getenv("HOME")+'/bin/'):
        default = yaml.load(open(os.getenv("HOME")+'/bin/default.aml'))
        # print "# Using FASTA file from ~/bin/default.aml"
        return default['FASTA_PATH']
    else:
        if 'FASTA_PATH' in os.environ:
            # print "# Using FASTA file from $FASTA_PATH variable"
            return os.environ['FASTA_PATH']
        else:
            exit('Provide FASTA file path using -g or setup default.aml in your bin folder')

def list_paths_in_current_dir(suffix=str(), stdin=False):
    """
    :param: suffix  -   lists paths in current directory ending with suffix only
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