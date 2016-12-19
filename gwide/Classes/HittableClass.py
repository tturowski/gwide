#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "1.0"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

import os, re, math, sys
import pandas as pd
import gwide.methods as gtk
import matplotlib.pyplot as plt

class HittableClass():
    def __init__(self, gtf, whole_name, n_rpM, out_prefix, read_stdin):
        self.gtf = gtf
        self.whole_name = whole_name
        self.n_rpM = n_rpM
        if out_prefix:
            self.out_prefix = out_prefix+'_'
        else:
            self.out_prefix = str()
        self.read_stdin = read_stdin

    def correlation(self, output, gene_class):
        print "# Calculate correlation is running..."

        no_of_reads = dict()
        genes_name = list()
        normalizator = float()
        classes = list()

        paths = gtk.list_paths_in_current_dir('hittable_reads.txt', stdin=self.read_stdin) #get paths of hittables
        experiments, paths = gtk.define_experiments(paths_in=paths, whole_name=self.whole_name) #extract experiments from paths

        data = pd.DataFrame(columns=[['gene_name', 'gene_id', 'type']+experiments]) # initialize Pandas DataFrame
        print "# Currently reading GTF file: "+ self.gtf
        #reading gtf line by line
        for line in open(self.gtf, 'r'):
            if not line.startswith('#'):
                line_elements = line.strip().split('\t')
                type = str(line_elements[1])
                if type not in classes:
                    classes.append(type)
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

        #create file with no of reads
        no_of_reads_file = open(self.out_prefix+"no_of_reads.table",'w')
        no_of_reads_file.write("# experiment"+'\t'+"mapped_reads"+'\t'+"total_reads"+'\n')

        #filling dataframe and close no_of_reads_file
        for path_no, path in enumerate(paths):
            name = experiments[path_no]
            no_of_reads[name] = dict()
            print "# Currently reading: "+path+" ..."
            for line in open(path, 'r'):
                if line.startswith('# total number of reads') and not line.startswith('# total number of reads without'):
                    total_reads = int(filter(str.isdigit, line)) # no of reads
                    no_of_reads[name]['total_reads'] = total_reads
                if line.startswith('# total number of single reads'):
                    total_mapped_reads = int(filter(str.isdigit, line)) # no of mapped reads
                    no_of_reads[name]['total_mapped_reads'] = total_mapped_reads
                    normalizator = 1000000.0/total_mapped_reads
                if not line.startswith('#'):
                    line_elements = line.strip().split('\t')
                    if len(line_elements) == 4: #needs modification if --rpkm applied for pyReadCounter
                        gene_name, hits = line_elements[0], int(line_elements[1])
                        if self.n_rpM == True:
                            data.loc[gene_name, name] = int(math.ceil(float(hits*normalizator)))
                        else:
                            data.loc[gene_name, name] = hits
            no_of_reads_file.write(name+'\t'+str(no_of_reads[name]['total_mapped_reads'])+'\t'+str(no_of_reads[name]['total_reads'])+'\n')
        no_of_reads_file.close()

        corr_dict = {"p" : "pearson" , "k" : "kendall" , "s" : "spearman"}
        if output == 'a':
            print "# Calculating all correlations..."
            for i in corr_dict:
                print "# Calculating correlations("+corr_dict[i]+")..."
                matrix = data.corr(method=corr_dict[i],min_periods=1)
                matrix.to_csv(self.out_prefix+"genome_wide_correlation_"+corr_dict[i]+".table", sep='\t')

        #calculate Pearson for different types
                if gene_class == True:
                    for this_type in classes:
                        new_data = data[data.type == this_type]
                        matrix = new_data.corr(method=corr_dict[i],min_periods=1)
                        matrix.to_csv(self.out_prefix+this_type+"_correlation_"+corr_dict[i]+".table", sep='\t')
        else:
            print "# Calculating correlations("+corr_dict[output]+")..."
            matrix = data.corr(method=corr_dict[output],min_periods=1)
            matrix.to_csv(self.out_prefix+"genome_wide_correlation_"+corr_dict[output]+".table", sep='\t')

            #calculate Pearson for different types
            if gene_class == True:
                for this_type in classes:
                    new_data = data[data.type == this_type]
                    matrix = new_data.corr(method=corr_dict[output],min_periods=1)
                    matrix.to_csv(self.out_prefix+this_type+"_correlation_"+corr_dict[output]+".table", sep='\t')

    def count(self, normalize=True, use_RPKM=False):
        no_of_reads = dict()
        genes_name = list()
        normalizator = float()

        paths = gtk.list_paths_in_current_dir('hittable_reads.txt', stdin=self.read_stdin) #get paths of hittables
        experiments, paths = gtk.define_experiments(paths_in=paths, whole_name=self.whole_name) #extract experiments from paths
        data = pd.DataFrame(columns=[['gene_name', 'gene_id', 'type']+experiments]) # initialize Pandas DataFrame

        #reading gtf file
        print "Currently reading GTF file: "+self.gtf
        for line in open(self.gtf, 'r'):
            if not line.startswith('#'):
                line_elements = line.strip().split('\t')
                type = str(line_elements[1])
                try:
                    gene_name = re.search("gene_name\s\"(.*?)\"", str(line_elements[8])).group(1)
                except:
                    gene_name = re.search("gene_id\s\"(.*?)\"", str(line_elements[8])).group(1) # when there is no gene name
                    print "No gene name in GTF file! Used gene id: "+gene_name+" as gene name."
                gene_id = re.search("gene_id\s\"(.*?)\"", str(line_elements[8])).group(1)
                if gene_name not in genes_name:
                    genes_name.append(gene_name)
                    gene_data = pd.DataFrame([[gene_name, gene_id, type]+([0]*len(experiments))], columns=(['gene_name', 'gene_id', 'type']+experiments))
                    data = data.append(gene_data, ignore_index=True)
        data = data.set_index(['gene_name'])

        #filling dataframe
        for path_no, path in enumerate(paths):
            name = experiments[path_no]
            print "Currently reading: "+path+"..."
            for line in open(path, 'r'):
                if line.startswith('# total number of single reads'):
                    total_mapped_reads = int(filter(str.isdigit, line)) # no of mapped reads
                    no_of_reads[name] = total_mapped_reads
                    if normalize == True: normalizator = 1000000.0/total_mapped_reads
                    else: normalizator = 1.0
                if not line.startswith('#'):
                    line_elements = line.strip().split('\t')
                    if len(line_elements) == 4:
                        gene_name, hits = line_elements[0], float(line_elements[1])
                        # print gene_name
                        data.loc[gene_name, name] = float(math.ceil(float(hits*normalizator)))
                    elif len(line_elements) == 6:
                        gene_name, hits, RPKM = line_elements[0], float(line_elements[1]), float(line_elements[2])
                        # print gene_name
                        if use_RPKM == False: data.loc[gene_name, name] = float(math.ceil(float(hits * normalizator)))
                        else: data.loc[gene_name, name] = float(math.ceil(float(RPKM * normalizator)))
        print "Creating output.tab file..."
        data.to_csv(self.out_prefix+'output.tab', sep='\t')

    def plot(self, print_single):

        paths = gtk.list_paths_in_current_dir('hittable_reads.txt', stdin=self.read_stdin) #get paths of hittables
        experiments, paths = gtk.define_experiments(paths_in=paths, whole_name=self.whole_name) #extract experiments from paths

        #initiating DataFrame
        data = pd.DataFrame(columns=[['group']+['legend']+experiments]) # initialize Pandas DataFrame
        data = data.set_index(['group'])
        general = dict()

        #filling DataFrame
        for path_no, hittable in enumerate(paths):
            name = experiments[path_no]
            general[name] = list() ## [total_mapped_reads, total_reads]
            for line in open(hittable):
                if line.startswith('# total number of mapped reads:'):
                    line_elements = line.strip().split('\t')
                    total_mapped_reads = int(line_elements[1])
                    general[name].append(total_mapped_reads)
                if line.startswith('# total number of reads') and not line.startswith('# total number of reads without'):
                    line_elements = line.strip().split('\t')
                    total_reads = int(line_elements[1])
                    general[name].append(total_reads)
                if line.startswith('##'):
                    line_elements = line.strip().split('\t')
                    type_of_reads = str(line_elements[0].strip('#').strip())
                    no_of_reads = int(line_elements[1])
                    data.loc[type_of_reads, name] = no_of_reads
            data = data.fillna(0)
        print data

        colors = ['lightblue', 'yellowgreen', 'darkred', 'gold',
                  'white','lightcoral','blue','pink', 'darkgreen',
                  'yellow','grey','violet','magenta','cyan']

        if print_single == False:
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            fig_no = 1
            plot_no = 1
            fig.add_subplot(3, 3, plot_no)
            plt.title('Legend:')
            plt.pie(data[experiments[0]], colors=colors, autopct='%1.1f%%', labeldistance=1.1, startangle=90)
            plt.legend(data.index, loc=0)
            for e in experiments:
                plot_no += 1
                fig.add_subplot(3, 3, plot_no)
                plt.pie(data[e], colors=colors, autopct='%1.1f%%', labeldistance=1.1, startangle=90)
                plt.axis('equal')
                plt.tight_layout()
                plt.title(e)
                plt.text(0,-0.9,'mapped reads/total reads: \n'+str(general[e][1])+'/'+str(general[e][0]),fontsize=12, horizontalalignment='center')

                if plot_no == 9:
                    plt.savefig('piecharts_'+str(fig_no)+'.png')
                    fig_no += 1
                    plt.clf()
                    plot_no = 0
            if plot_no > 0:
                plt.savefig('piecharts_'+str(fig_no)+'.png')
                plt.clf()

        elif print_single == True:
            print "Plotting piecharts in separate files..."
            for e in experiments:
                plt.pie(data[e], colors=colors, autopct='%1.1f%%', labeldistance=1.1, startangle=90)
                plt.legend(data.index, loc=4)
                plt.axis('equal')
                plt.tight_layout()
                plt.title(e)
                plt.text(-0.9,-0.9,'mapped reads/total reads: \n'+str(general[e][1])+'/'+str(general[e][0]),fontsize=12, horizontalalignment='center')
                plt.savefig(e+'.png')
                plt.clf()
        print 'Done.'