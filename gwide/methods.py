import os, yaml, sys, re
import numpy as np, numpy.random
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
        file_path = path.split('/')
        file_name = file_path[len(file_path)-1]
        if whole_name == False: name = "_".join(file_name.split('_')[0:3]) #take fist three elements of file name as experiment name
        elif whole_name == True: name = file_name.rstrip('_hittable_reads.txt')
        experiments.append(name)
        if len(experiments) != len(paths):
            exit("No. of experiments is not equal to no. of paths")
    return experiments, paths

def calGC(dataset=pd.DataFrame(), calFor=['G','C']):
    '''
    Returns GC content in a given dataset
    :param dataset: Pandas DataFrame with "nucleotide" column
    :return: fraction of GC content
    '''
    return float(len(dataset[dataset['nucleotide'].isin(calFor)]))/float(len(dataset))


def expNameParser(name, additional_tags=list(), order='b_d_e_p'):
    '''
    Function handles experiment name; recognizes AB123456 as experiment date; BY4741 or HTP or given string as bait protein
    :param name:
    :param additional_tags: list of tags
    :param output: default 'root' ; print other elements when 'all'
    :param order: defoult 'b_d_e_p' b-bait; d-details, e-experiment, p-prefix
    :return: reordered name
    '''
    tag_list = ['HTP', 'HTG', 'HTF', 'BY4741'] + additional_tags
    output_dict = {'b': str(), 'd': str(), 'e': str(), 'p': list()}  # bait; details; experiment; prefix
    name_elements = name.split('_')
    for e in name_elements:
        tag_in_e = [tag for tag in tag_list if tag in e]
        if tag_in_e and len(tag_in_e) >= 1:
            output_dict['b'] = e  # bait
            try:
                output_dict['d'] = name.split(tag_in_e[0], 1)[1].strip('_')  # details
            except:
                output_dict['d'] = 'wt'
                print 'WARNING: wt added for '+name
        elif re.search(r"[a-zA-Z][a-zA-Z]\d{6}", e) or re.search(r"[a-zA-Z][a-zA-Z][a-zA-Z]\d{6}", e):
            output_dict['e'] = e  # experiment name
            try:
                output_dict['p'] = name.split(e, 1)[0].strip('_')  # prefix
            except:
                output_dict['p'] = ''

    if len(output_dict['b']) < 3 or len(output_dict['e']) < 8:
        print output_dict
        sys.exit("ERROR: Can not define experiment or bait for "+name)

    return_list = list()
    for out in order.split('_'):
        return_list.append(output_dict[out])

    return '_'.join(return_list).strip('_')

def cleanNames(df=pd.DataFrame(), additional_tags=list()):
    '''Cleans some problems with names if exist'''
    for tag in additional_tags: df.columns = [f.replace(tag, tag+'HTP') for f in list(df.columns.values)]
    df.columns = [f.replace('HTPHTP', 'HTP').replace('HTPHTP', 'HTP') for f in list(df.columns.values)]
    return df

def indexOrder(df=pd.DataFrame(), additional_tags=list(), output='root', order='b_d_e_p'):
    '''Aplly expNameParser to whole dataframe
    :param order: defoult 'b_d_e_p' b-bait; d-details, e-experiment, p-prefix
    '''
    df = cleanNames(df, additional_tags=additional_tags)
    df.columns = [expNameParser(f, additional_tags=additional_tags, order=order) for f in list(df.columns.values)]
    return df.sort_index(axis=1)

def filterExp(datasets, let_in=[''], let_out=['wont_find_this_string']):
    '''for pd.DataFrame() or dict(). Returns object with filtered columns/keys.

    datasets : DataFrame() or dict()
      DataFrame() or dict() with exp name as a key

    let_in : list()
      list() with elements of name to filter in

    let_out : list()
      list() with elements of name to filter out

    Returns
    -------
    DataFrame() or dict()
    '''
    #for pd.DataFrame()
    if type(datasets) == type(pd.DataFrame()):
        output_df = pd.DataFrame()
        for f in [d for d in list(datasets.columns.values) if all(i in d for i in let_in) and all(o not in d for o in let_out)]:
            output_df[f]=datasets[f]
        return output_df
    #for dict()
    elif type(datasets) == type(dict()):
        output_dict = dict()
        for f in [d for d in list(datasets.keys()) if all(i in d for i in let_in) and all(o not in d for o in let_out)]:
            output_dict[f]=datasets[f]
        return output_dict


def calculateFDR(data=pd.Series(), iterations=100, target_FDR=0.05):
    '''calculates False Discovery Rate (FDR) for a given dataset.

    data : pd.Series()

    iterations : int()
        number of iterations. Default = 100
    target_FDR : float()
        Detault = 0.05
    Returns
    -----------
    Series()

    '''
    normalized_data = data / data.sum()  # normalize data
    random_df = pd.DataFrame(np.random.rand(len(data), iterations))  # generating random datasets
    normalized_random_df = random_df / random_df.sum()  # normalize random df

    # getting Discovery Rate (DR)
    DR_df = normalized_random_df.T >= normalized_data  # random discovery rate
    DR = DR_df.sum() / iterations

    # comparing DR to False Discovery Rate (FDR)
    FDR = DR <= target_FDR
    return data * FDR

