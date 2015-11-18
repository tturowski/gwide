__author__ = 'tturowski'

import os, yaml
import numpy as np
import pandas as pd

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
