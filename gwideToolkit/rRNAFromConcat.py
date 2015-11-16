#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "1.3"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

import numpy as np
import sys, collections, re, os
from pypeaks import Data
from pyCRAC.Parsers import GTF2
import matplotlib.pyplot as plt
import pandas as pd

class rRNAFromConcat():
    def __init__(self, gtf_file, five_prime_flank, three_prime_flank, hits_threshold, lookahead, prefix, print_peaks, normalized):
        self.gtf_file           =   str(gtf_file)
        self.gtf                =   GTF2.Parse_GTF()
        self.gtf.read_GTF(self.gtf_file)
        self.genes              =   dict()
        self.data               =   dict()
        self.id_to_names        =   dict()
        self.genes_name_list    =   list()
        self.genes_id_list      =   list()
        self.five_prime_flank   =   five_prime_flank
        self.three_prime_flank  =   three_prime_flank
        self.hits_threshold     =   hits_threshold
        self.lookahead          =   lookahead
        self.prefix             =   str(prefix)
        self.print_peaks        =   print_peaks
        self.experiments        =   list()
        self.normalized         =   normalized  # -n option allows for using normalized dataset (reads p[er Million)

    def read_concat_file(self, concat_file, null_substitution):
        print '# Reading concat file.'
        rRNA_list = ['RDN37-1', 'RDN37-2']
        for line in concat_file:
            try:
                if line[0] == "#": continue
                line_elements = line.strip().split('\t')
                if self.normalized == False:
                    gene_name, position, nucleotide, hits, experiment = line_elements[0],int(line_elements[1]), line_elements[2], int(line_elements[3]), line_elements[6]
                elif self.normalized == True:
                    gene_name, position, nucleotide, hits, experiment = line_elements[0],int(line_elements[1]), line_elements[2], int(line_elements[7]), line_elements[6]
                if gene_name not in rRNA_list: continue
        #adding new entry in genes dict
                if gene_name not in self.genes:
                    gene_id = self.gtf.genes[gene_name]['gene_id']
                    gene_length = self.gtf.geneLength(gene_name)
                    self.genes[gene_name] = {
                        'gene_name'     :   gene_name,
                        'gene_id'       :   gene_id,
                        'gene_length'   :   gene_length,
                        }
                    self.id_to_names[gene_id] = gene_name
                    self.genes_name_list.append(gene_name)
                    self.genes_id_list.append(gene_id)
        #checking does experiment was noted before
                if gene_name not in self.data:
                    self.data[gene_name] = dict()
                if experiment not in self.data[gene_name]:
                    frame_index = range(1,(self.five_prime_flank + self.genes[gene_name]['gene_length'] + self.three_prime_flank+1))
                    self.data[gene_name][experiment] = pd.DataFrame(index=frame_index, columns=['position','nucleotides', 'exon', 'RT'])
                if experiment not in self.experiments:
                    self.experiments.append(experiment)
        #filling DataFrame
                # if position <= self.five_prime_flank :
                self.data[gene_name][experiment].loc[position, 'position'] = position - self.five_prime_flank + 700
                self.data[gene_name][experiment].loc[position, 'nucleotides'] = nucleotide
                if null_substitution == True and hits == 0: ## option necessry for ratio and log calculation
                    self.data[gene_name][experiment].loc[position, 'hits'] = 1
                else:
                    self.data[gene_name][experiment].loc[position, 'hits'] = hits

            except IndexError:
                sys.stderr.write("\nIndexError at line:\n%s\n" % line)
                pass

        self.experiments.sort()

    def find_peaks(self):
        print '# Finding peaks.'
        for i in self.data:
            for e in self.data[i]:
                hist = Data(list(self.data[i][e].index), list(self.data[i][e]['hits']), smoothness=1, default_smooth=False)
                hist.normalize()
                try:
                    hist.get_peaks(method="slope", peak_amp_thresh=0.00005, valley_thresh=0.00003, intervals=None,
                                   lookahead=self.lookahead, avg_interval=100)
                    self.data[i][e]['peaks'] = sorted(np.array(hist.peaks['peaks'][0]).tolist())
                    self.data[i][e]['valleys'] = sorted(np.array(hist.peaks['valleys'][0]).tolist())
                except ValueError:
                    self.data[i][e]['peaks'] = []
                    self.data[i][e]['valleys'] = []
        return True

# ### preparing dataframes to for plotting
#     def create_features_to_plot(self):
#         for i in self.data:
#             frame_index = range(1,(self.five_prime_flank + self.genes[i]['gene_length'] + self.three_prime_flank+1))
#
#             for e in self.data[i]:
#                 histogram = list(self.data[i][e]['hits'])
#                 hist_max = max(histogram)
#     # histogram for peaks
#                 if self.print_peaks == True:
#                     peaks = self.data[i][e]['peaks']
#                     for c in frame_index:
#                         if c in peaks:
#                             self.data[i][e].loc[c, 'peaks_to_plot'] = hist_max
#                         else:
#                             self.data[i][e].loc[c, 'peaks_to_plot'] = 0
#
#     def slice_dataframe(self):
#         five_prime_to_out = self.five_prime_flank - 50 # maybe changed - it simply gives only -50 flank independly what -r value was used for pyPileups
#         for i in self.data:
#             for e in self.data[i]:
#                 self.data[i][e] = self.data[i][e][five_prime_to_out::]
#         return True

    # get smaller array -> -600 on the 5` end
    def slice_data(self):
        for i in self.data:
            for e in self.data[i]:
                self.data[i][e] = self.data[i][e][:-600:]

### making output plots, one gene, different experiments per page
    def print_rRNA(self):
        print '# Plotting RDN37 (all experiments). Did you set up both flanks to 1000 bp?'
        fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
        fig_no = 0
        for i_gene_id in self.genes_id_list:
            gene_name = self.id_to_names[i_gene_id]
            plot_no = 0
            for e in self.experiments:
                plot_no += 1
                fig.add_subplot(3, 1, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                if self.print_peaks == True:
                    plt.plot(self.data[gene_name][e]['peaks_to_plot'])

                first_flank = 300
                ETS = 700
                eighteen = 1800
                its_one = 361
                five_dot_eight = 158
                its_two = 232
                twenty_five = 3396
                ETS_end = 211

                sliced = collections.OrderedDict()
                sliced['5flank'] = self.data[gene_name][e][0:first_flank:]       #5` flank
                sliced['5ETS'] = self.data[gene_name][e][first_flank:(first_flank+ETS):]    #5` ETS
                sliced['18S'] = self.data[gene_name][e][(first_flank+ETS):(first_flank+ETS+eighteen):]   #18S
                sliced['ITS1'] = self.data[gene_name][e][(first_flank+ETS+eighteen):(first_flank+ETS+eighteen+its_one):]   #ITS1
                sliced['5.8S'] = self.data[gene_name][e][(first_flank+ETS+eighteen+its_one):(first_flank+ETS+eighteen+its_one+five_dot_eight):]   #5.8S
                sliced['ITS2'] = self.data[gene_name][e][(first_flank+ETS+eighteen+its_one+five_dot_eight):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):]   #ITS2
                sliced['25S'] = self.data[gene_name][e][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):]   #25S
                sliced['3ETS`'] = self.data[gene_name][e][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end):]   #3` ETS
                sliced['3flank'] = self.data[gene_name][e][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end)::]       #3` flank
                for part in sliced:
                    x_array = np.array(sliced[part]['position'])
                    y_array = np.array(sliced[part]['hits'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array, label=part)
                    if part in ['18S', '5.8S', '25S']:
                        plt.axvspan(min(sliced[part]['position']), max(sliced[part]['position']), alpha=0.2, color='#0892d0')
                plt.grid()
                legend = plt.legend(loc='upper right', shadow=True, fontsize=10)

                plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                if plot_no == 3:
                    plot_no = 0
                    fig_no += 1
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png')
                    plt.clf()
            if plot_no > 0:
                plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    ###function to plot ration between sample with -a parameter and -b parameter
    def fig_ratio(self, a, b):
        print '# Plotting ratio for '+a+' divided by '+b+' (all experiments).'
        new_exp_list = self.group_experiments(a,b)
        fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
        fig_no = 0
        for i_gene_id in self.genes_id_list:
            gene_name = self.id_to_names[i_gene_id]
            plot_no = 0
            for e in new_exp_list:
                plot_no += 1
                fig.add_subplot(3, 1, plot_no)
                plt.tight_layout()
                plt.title(e[0]+'/'+e[1]+' ratio')
                plt.ylabel("no. of reads")
                self.data[gene_name][e[0]]['ratio'] = self.data[gene_name][e[0]]['hits']/self.data[gene_name][e[1]]['hits'] ###getting ratio

                first_flank = 300
                ETS = 700
                eighteen = 1800
                its_one = 361
                five_dot_eight = 158
                its_two = 232
                twenty_five = 3396
                ETS_end = 211

                sliced = collections.OrderedDict()
                sliced['5flank'] = self.data[gene_name][e[0]][0:first_flank:]       #5` flank
                sliced['5ETS'] = self.data[gene_name][e[0]][first_flank:(first_flank+ETS):]    #5` ETS
                sliced['18S'] = self.data[gene_name][e[0]][(first_flank+ETS):(first_flank+ETS+eighteen):]   #18S
                sliced['ITS1'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen):(first_flank+ETS+eighteen+its_one):]   #ITS1
                sliced['5.8S'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one):(first_flank+ETS+eighteen+its_one+five_dot_eight):]   #5.8S
                sliced['ITS2'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):]   #ITS2
                sliced['25S'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):]   #25S
                sliced['3ETS`'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end):]   #3` ETS
                sliced['3flank'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end)::]       #3` flank
                for part in sliced:
                    x_array = np.array(sliced[part]['position'])
                    y_array = np.array(sliced[part]['ratio'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array, label=part)
                    if part in ['18S', '5.8S', '25S']:
                        plt.axvspan(min(sliced[part]['position']), max(sliced[part]['position']), alpha=0.2, color='#0892d0')
                plt.grid()
                legend = plt.legend(loc='upper right', shadow=True, fontsize=10)

                plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                if plot_no == 3:
                    plot_no = 0
                    fig_no += 1
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_ratio_'+str(fig_no)+'.png')
                    plt.clf()
            if plot_no > 0:
                plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_ratio_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    ###function to plot ration between sample with -a parameter and -b parameter
    def fig_log2ratio(self, a, b):
        print '# Plotting log2 ratio for '+a+' divided by '+b+' (all experiments).'
        new_exp_list = self.group_experiments(a,b)
        fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
        fig_no = 0
        for i_gene_id in self.genes_id_list:
            gene_name = self.id_to_names[i_gene_id]
            plot_no = 0
            for e in new_exp_list:
                plot_no += 1
                fig.add_subplot(3, 1, plot_no)
                plt.tight_layout()
                plt.title(e[0]+'/'+e[1]+' log2 ratio')
                plt.ylabel("no. of reads")
                self.data[gene_name][e[0]]['log2'] = np.log2(self.data[gene_name][e[0]]['hits']/self.data[gene_name][e[1]]['hits']) ###getting ratio
                self.data[gene_name][e[0]].to_csv(sep='\t', path_or_buf=(os.getcwd()+'/'+e[0]+'.tab'))
                self.data[gene_name][e[0]].to_csv(sep=',', path_or_buf=(os.getcwd()+'/'+e[0]+'.csv'))

                first_flank = 300
                ETS = 700
                eighteen = 1800
                its_one = 361
                five_dot_eight = 158
                its_two = 232
                twenty_five = 3396
                ETS_end = 211

                sliced = collections.OrderedDict()
                sliced['5flank'] = self.data[gene_name][e[0]][0:first_flank:]       #5` flank
                sliced['5ETS'] = self.data[gene_name][e[0]][first_flank:(first_flank+ETS):]    #5` ETS
                sliced['18S'] = self.data[gene_name][e[0]][(first_flank+ETS):(first_flank+ETS+eighteen):]   #18S
                sliced['ITS1'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen):(first_flank+ETS+eighteen+its_one):]   #ITS1
                sliced['5.8S'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one):(first_flank+ETS+eighteen+its_one+five_dot_eight):]   #5.8S
                sliced['ITS2'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):]   #ITS2
                sliced['25S'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):]   #25S
                sliced['3ETS`'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end):]   #3` ETS
                sliced['3flank'] = self.data[gene_name][e[0]][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end)::]       #3` flank
                for part in sliced:
                    x_array = np.array(sliced[part]['position'])
                    y_array = np.array(sliced[part]['log2'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array, label=part)
                    if part in ['18S', '5.8S', '25S']:
                        plt.axvspan(min(sliced[part]['position']), max(sliced[part]['position']), alpha=0.2, color='#0892d0')
                plt.grid()
                legend = plt.legend(loc='upper right', shadow=True, fontsize=10)

                plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                if plot_no == 3:
                    plot_no = 0
                    fig_no += 1
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_log2ratio_'+str(fig_no)+'.png')
                    plt.clf()
            if plot_no > 0:
                plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_log2ratio_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    def group_experiments(self, a, b):
        print 'Looking for groups of experiments:'
### function identifiying 'to_divide' and 'divisor' for ratio results plotting.
### Please use nomenclature gene_variant
### IMPORTANT: only two experiments may contain the same root(gene) added to 'a' or 'b'
### i.e. list_of_experiments = ['A190_MD', 'A190_total', 'A135_MD', 'A135_total'] additional A190 or A135 is forbidden
        a_experiments = list()
        b_experiments = list()
        for e in self.experiments: ##making list with experiments containing -a or -b
            if re.search(a, e):
                a_experiments.append(e)
            elif re.search(b, e):
                b_experiments.append(e)
            else:
                print "fig_ratio module found experiment "+e+" neither containing -a nor -b parameter. \n " \
                                                             "This experiment won't be considered"

        if len(a_experiments) == len(b_experiments): ##checking no. of experiments to divide each other
            pass
        else:
            print "No of experiments with -a and -b parameter is unequal:"
            print "-a exp: "+','.join(a_experiments)+"\n -b exp: "+','.join(b_experiments)
            exit()

        a_experiments.sort()
        b_experiments.sort()
        ratio_exp_list = list()
        for e in a_experiments:
            gene = re.search(r'(.*)_'+a, e).group(1) ### get prefix with name of gene
            second = [d for d in b_experiments if re.search(r'^'+gene,d)]
            if len(second) > 1:
                print 'To many elements to divide for '+gene
                exit()
            couple = [e, second[0]]
            ratio_exp_list.append(couple)
        print 'Found pairs:'
        print ratio_exp_list
        return ratio_exp_list
####################

    def test_print(self, what):
        pd.set_option('display.max_rows', 5000)
        print what
        pd.reset_option('display.max_rows')
        return True