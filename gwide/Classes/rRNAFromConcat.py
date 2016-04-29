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
        self.npM                =   normalized  # -n option allows for using normalized dataset (reads p[er Million)

    # def read_concat_file(self, concat_file, null_substitution):
    #     print '# Reading concat file.'
    #     rRNA_list = ['RDN37-1', 'RDN37-2']
    #     for line in concat_file:
    #         try:
    #             if line[0] == "#": continue
    #             line_elements = line.strip().split('\t')
    #             if self.normalized == False:
    #                 gene_name, position, nucleotide, hits, experiment = line_elements[0],int(line_elements[1]), line_elements[2], int(line_elements[3]), line_elements[6]
    #             elif self.normalized == True:
    #                 gene_name, position, nucleotide, hits, experiment = line_elements[0],int(line_elements[1]), line_elements[2], int(line_elements[7]), line_elements[6]
    #             if gene_name not in rRNA_list: continue
    #     #adding new entry in genes dict
    #             if gene_name not in self.genes:
    #                 gene_id = self.gtf.genes[gene_name]['gene_id']
    #                 gene_length = self.gtf.geneLength(gene_name)
    #                 self.genes[gene_name] = {
    #                     'gene_name'     :   gene_name,
    #                     'gene_id'       :   gene_id,
    #                     'gene_length'   :   gene_length,
    #                     }
    #                 self.id_to_names[gene_id] = gene_name
    #                 self.genes_name_list.append(gene_name)
    #                 self.genes_id_list.append(gene_id)
    #     #checking does experiment was noted before
    #             if gene_name not in self.data:
    #                 self.data[gene_name] = dict()
    #             if experiment not in self.data[gene_name]:
    #                 frame_index = range(1,(self.five_prime_flank + self.genes[gene_name]['gene_length'] + self.three_prime_flank+1))
    #                 self.data[gene_name][experiment] = pd.DataFrame(index=frame_index, columns=['position','nucleotides', 'exon', 'RT'])
    #             if experiment not in self.experiments:
    #                 self.experiments.append(experiment)
    #     #filling DataFrame
    #             # if position <= self.five_prime_flank :
    #             self.data[gene_name][experiment].loc[position, 'position'] = position - self.five_prime_flank + 700
    #             self.data[gene_name][experiment].loc[position, 'nucleotides'] = nucleotide
    #             if null_substitution == True and hits == 0: ## option necessry for ratio and log calculation
    #                 self.data[gene_name][experiment].loc[position, 'hits'] = 1
    #             else:
    #                 self.data[gene_name][experiment].loc[position, 'hits'] = hits
    #         except IndexError:
    #             sys.stderr.write("\nIndexError at line:\n%s\n" % line)
    #             pass
    #     self.experiments.sort()

    def read_csv(self, concat_file):
        print "# Reading CSV file..."
        header = ['gene','position','nucleotide','hits','substitutions','deletions','exp_name','n_hits','n_substitutions','n_deletions']
        concat_csv = pd.read_csv(concat_file, sep='\t', names=header, comment='#')
        self.experiments = sorted(set(concat_csv['exp_name'].tolist()))
        self.genes_name_list = sorted(set(concat_csv['gene'].tolist()))
        #adding new entry in genes dict
        print "# Getting details from GTF file and filling dataframe..."
        for gene_name in self.genes_name_list:
            gene_id = self.gtf.genes[gene_name]['gene_id']
            gene_length = self.gtf.geneLength(gene_name)
            self.genes[gene_name] = {
                'gene_name'     :   gene_name,
                'gene_id'       :   gene_id,
                'gene_length'   :   gene_length,
                                    }
            self.id_to_names[gene_id] = gene_name
            self.genes_id_list.append(gene_id)
        #filling dataframe
            gene_csv = concat_csv[concat_csv.gene == gene_name].set_index(['position'])
            frame_index = range(1,(self.five_prime_flank + gene_length + self.three_prime_flank+1))
            positions = range(-self.five_prime_flank,0)+range(1,(gene_length+self.three_prime_flank+1))
            # columns = ['position','nucleotide']+[self.experiments]
            columns = ['position','nucleotide']
            self.data[gene_name] = pd.DataFrame(index=frame_index, columns=columns)
            self.data[gene_name]['position'] = positions
            self.data[gene_name]['nucleotide'] = gene_csv['nucleotide'][:(self.five_prime_flank + gene_length + self.three_prime_flank):]
            for e in self.experiments:
                self.data[gene_name][e] = gene_csv[gene_csv.exp_name == e]['hits']
                if self.npM == True:
                    self.data[gene_name][(e+"_nrpm")] = gene_csv[gene_csv.exp_name == e]['n_hits']    #to work with data normalized reads per M
            self.data[gene_name] = self.data[gene_name].fillna(0)
            # print self.data[gene_name][230:330]
        self.genes_id_list.sort()
        return True

    def calculate(self, details=False, ntotal=False, nmax=False, pscounts=False):
        if details == True or nmax == True or ntotal == True:
            print '# Calculating readthrough and others...'
        elif details == True and len(self.experiments) > 1:
            print '# WORNING: -d parameter works only with one experiment.'
            exit()
        else:
            print '# Calculating readthrough...'

        for gene_name in self.genes:
        #     transcription_start =   self.five_prime_flank-21
        #     gene_length         =   self.genes[gene_name]['gene_length']
        #     gene_end            =   self.five_prime_flank + gene_length
        #     RT_begin            =   self.five_prime_flank + gene_length + self.readthrough_start
        #     gene_middle         =   self.five_prime_flank + ( gene_length / 2 )
        #     three_middle           =   self.five_prime_flank + gene_length + (self.three_prime_flank / 2)
        #     three_one_third        =   self.five_prime_flank + gene_length + (self.three_prime_flank / 3)
        #     three_two_third        =   self.five_prime_flank + gene_length + (2*(self.three_prime_flank / 3))
        #
        # #getting intron length
        #     introns = list()
        #     if not self.genes[gene_name]['introns'][0]:
        #         intron_length = 0
        #     else:
        #         for intron in range(0,len(self.genes[gene_name]['introns'][0])):
        #             introns.append(str(self.genes[gene_name]['introns'][0][intron]))
        #         intron_length = int(''.join(map(str,introns)))
        #         intron_start_stop = self.genes[gene_name]['introns'][1][0]      #start and stop of first intron only!
            for exp in self.experiments:
        # !! changes in exp name !!
                exp_old = exp
                # try:
                if max(list(self.data[gene_name][exp])) >= self.hits_threshold:
            #normalization options
                    if nmax == True:
                        exp = exp_old+'_nmax'
                        self.data[gene_name][exp] = self.data[gene_name][exp_old]/self.data[gene_name][exp_old].max()
                        if pscounts == True:
                            self.data[gene_name][exp] = self.data[gene_name][exp].add(0.000001) #adding pseudocounts
                    if ntotal == True:
                        exp = exp_old+'_ntotal'
                        self.data[gene_name][exp] = self.data[gene_name][exp_old]/self.data[gene_name][exp_old].sum()
                        if pscounts == True:
                            self.data[gene_name][exp] = self.data[gene_name][exp].add(0.000001) #adding pseudocounts
                    if pscounts == True:
                        self.data[gene_name][exp_old] = self.data[gene_name][exp_old].add(10) #adding pseudocounts

                # #slicing dataframes
                #         total       =   self.data[gene_name][transcription_start:].sum()[exp]
                #         total_av    =   self.data[gene_name][transcription_start:].mean()[exp]
                #         total_med   =   self.data[gene_name][transcription_start:].median()[exp]
                #         total_SD    =   self.data[gene_name][transcription_start:].std()[exp]
                #         c           =   float(self.data[gene_name][RT_begin:].sum()[exp])
                #         if details==True:
                #             g       =   self.data[gene_name][transcription_start:gene_end].sum()[exp]
                #             g_av    =   self.data[gene_name][transcription_start:gene_end].mean()[exp]
                #             g_med   =   self.data[gene_name][transcription_start:gene_end].median()[exp]
                #             g_SD    =   self.data[gene_name][transcription_start:gene_end].std()[exp]
                #             a1      =   float(self.data[gene_name][transcription_start:gene_middle].sum()[exp])
                #             a2      =   self.data[gene_name][gene_middle+1:gene_end].sum()[exp]
                #             d       =   self.data[gene_name][gene_end+1:].sum()[exp]
                #             d_av    =   self.data[gene_name][gene_end+1:].mean()[exp]
                #             d_med   =   self.data[gene_name][gene_end+1:].median()[exp]
                #             d_SD    =   self.data[gene_name][gene_end+1:].std()[exp]
                #             e1      =   float(self.data[gene_name][gene_end+1:three_middle].sum()[exp])
                #             e2      =   self.data[gene_name][three_middle+1:].sum()[exp]
                #             f1      =   float(self.data[gene_name][gene_end+1:three_one_third].sum()[exp])
                #             f2      =   self.data[gene_name][three_one_third+1:three_two_third].sum()[exp]
                #             f3      =   self.data[gene_name][three_two_third+1:].sum()[exp]
                #             if intron_length > 0:
                #                 b1  =   float(self.data[gene_name][transcription_start:intron_start_stop[0]+self.five_prime_flank].sum()[exp])
                #                 b2  =   float(self.data[gene_name][intron_start_stop[0]+1+self.five_prime_flank:intron_start_stop[1]+self.five_prime_flank].sum()[exp])
                #                 b3  =   self.data[gene_name][intron_start_stop[1]+1+self.five_prime_flank:gene_end].sum()[exp]
                # #calculating
                #         RT = np.float64(c) / total #allows for dividing by 0
                #         if details==True:
                #             a = np.float64(a1) / a2
                #             e = np.float64(e1) / e2
                #             f = np.float64(f1) / f3
                #
                #             if intron_length > 0:
                #                 b = np.float64(b1) / b3
                #                 i = np.float64(b2) / (b1 + b3)
                #             else:
                #                 b = 0
                #                 i = 0
                #     else:
                #         RT, total, total_av, total_med, total_SD, c = ['too_low_reads'] * 6
                #         if details == True:
                #             g, g_av, g_med, g_SD, d, d_av, d_med, d_SD, a, e, f, b, i = ['too_low_reads'] * 13
                # except KeyError as e:
                    # print "Error raised by key: "+str(e)
                #     RT, total, total_av, total_med, total_SD, c = ['no_reads'] * 6
                #     if details == True:
                #             g, g_av, g_med, g_SD, d, d_av, d_med, d_SD, a, e, f, b, i = ['no_reads'] * 13
                # #stroing in dictionary
                # self.genes[gene_name]['RT'][exp_old] = RT
                # if details==True:
                #     self.genes[gene_name]['total']      =   total
                #     self.genes[gene_name]['total_av']   =   total_av
                #     self.genes[gene_name]['total_med']  =   total_med
                #     self.genes[gene_name]['total_std']  =   total_SD
                #     self.genes[gene_name]['a']          =   a
                #     self.genes[gene_name]['b']          =   b
                #     self.genes[gene_name]['i']          =   i
                #     self.genes[gene_name]['e']          =   e
                #     self.genes[gene_name]['f']          =   f
                #     self.genes[gene_name]['d']      =   d
                #     self.genes[gene_name]['d_av']   =   d_av
                #     self.genes[gene_name]['d_med']  =   d_med
                #     self.genes[gene_name]['d_std']  =   d_SD
                #     self.genes[gene_name]['g']      =   g
                #     self.genes[gene_name]['g_av']   =   g_av
                #     self.genes[gene_name]['g_med']  =   g_med
                #     self.genes[gene_name]['g_std']  =   g_SD
        return True

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
        fig = plt.figure(figsize=(12, 9), dpi=300, facecolor='w', edgecolor='k')
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
                sliced['5flank'] = self.data[gene_name][0:first_flank:]       #5` flank
                sliced['5ETS'] = self.data[gene_name][first_flank:(first_flank+ETS):]    #5` ETS
                sliced['18S'] = self.data[gene_name][(first_flank+ETS):(first_flank+ETS+eighteen):]   #18S
                sliced['ITS1'] = self.data[gene_name][(first_flank+ETS+eighteen):(first_flank+ETS+eighteen+its_one):]   #ITS1
                sliced['5.8S'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one):(first_flank+ETS+eighteen+its_one+five_dot_eight):]   #5.8S
                sliced['ITS2'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):]   #ITS2
                sliced['25S'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):]   #25S
                sliced['3ETS`'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end):]   #3` ETS
                sliced['3flank'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end)::]       #3` flank
                for part in sliced:
                    # print sliced[part]
                    x_array = np.array(sliced[part]['position'])
                    y_array = np.array(sliced[part][e])
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
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png', dpi=300)
                    plt.clf()
            if plot_no > 0:
                plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no+1)+'.png', dpi=300)
                plt.clf()
        return True

### making output plots, one gene, different experiments per page
    def single_rRNA(self):
        print '# Plotting RDN37 (all experiments). Did you set up both flanks to 1000 bp?'
        fig = plt.figure(figsize=(12, 9), dpi=300, facecolor='w', edgecolor='k')
        fig_no = 0
        for i_gene_id in self.genes_id_list:
            gene_name = self.id_to_names[i_gene_id]
            for e in self.experiments:
                fig.add_subplot(3, 1, 1)
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
                sliced['5flank'] = self.data[gene_name][0:first_flank:]       #5` flank
                sliced['5ETS'] = self.data[gene_name][first_flank:(first_flank+ETS):]    #5` ETS
                sliced['18S'] = self.data[gene_name][(first_flank+ETS):(first_flank+ETS+eighteen):]   #18S
                sliced['ITS1'] = self.data[gene_name][(first_flank+ETS+eighteen):(first_flank+ETS+eighteen+its_one):]   #ITS1
                sliced['5.8S'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one):(first_flank+ETS+eighteen+its_one+five_dot_eight):]   #5.8S
                sliced['ITS2'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):]   #ITS2
                sliced['25S'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):]   #25S
                sliced['3ETS`'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five):(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end):]   #3` ETS
                sliced['3flank'] = self.data[gene_name][(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end)::]       #3` flank
                for part in sliced:
                    # print sliced[part]
                    x_array = np.array(sliced[part]['position'])
                    y_array = np.array(sliced[part][e])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array, label=part)
                    if part in ['18S', '5.8S', '25S']:
                        plt.axvspan(min(sliced[part]['position']), max(sliced[part]['position']), alpha=0.2, color='#0892d0')
                plt.grid()
                legend = plt.legend(loc='upper right', shadow=True, fontsize=10)
                fig_no += 1
                plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png', dpi=300)
                plt.clf()
        return True

    ###function to plot ration between sample with -a parameter and -b parameter
    def fig_ratio(self, a, b):
        print '# Plotting ratio for '+a+' divided by '+b+' (all experiments).'
        new_exp_list = self.group_experiments(a,b)
        fig = plt.figure(figsize=(12, 9), dpi=300, facecolor='w', edgecolor='k')
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
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_ratio_'+str(fig_no)+'.png', dpi=300)
                    plt.clf()
            if plot_no > 0:
                plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_ratio_'+str(fig_no+1)+'.png', dpi=300)
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

    def correlations(self, output='a', elements=True, ntotal=True):
        print '# Calculating RDN37 plot correlations. Did you set up both flanks to 1000 bp?'
        for i_gene_id in self.genes_id_list:
            gene_name = self.id_to_names[i_gene_id]
            #drop additional collumns
            self.data[gene_name] = self.data[gene_name].fillna(0)
            self.data[gene_name] = self.data[gene_name].drop('position', axis=1)
            self.data[gene_name] = self.data[gene_name].drop('nucleotide', axis=1)

            first_flank = 300+1
            ETS = 700
            eighteen = 1800
            its_one = 361
            five_dot_eight = 158
            its_two = 232
            twenty_five = 3396
            ETS_end = 211

            self.data[gene_name]['element'] = 'unknown'
            # self.data[gene_name].update(pd.DataFrame({'element' : ['5flank']}, index=range(0,first_flank)))
            self.data[gene_name].update(pd.DataFrame({'element' : ['5ETS']}, index=range(first_flank,first_flank+ETS)))
            self.data[gene_name].update(pd.DataFrame({'element' : ['18S']}, index=range(first_flank+ETS,first_flank+ETS+eighteen)))
            self.data[gene_name].update(pd.DataFrame({'element' : ['ITS1']}, index=range(first_flank+ETS+eighteen,first_flank+ETS+eighteen+its_one)))
            self.data[gene_name].update(pd.DataFrame({'element' : ['5.8S']}, index=range(first_flank+ETS+eighteen+its_one,first_flank+ETS+eighteen+its_one+five_dot_eight)))
            self.data[gene_name].update(pd.DataFrame({'element' : ['ITS2']}, index=range(first_flank+ETS+eighteen+its_one+five_dot_eight,first_flank+ETS+eighteen+its_one+five_dot_eight+its_two)))
            self.data[gene_name].update(pd.DataFrame({'element' : ['25S']}, index=range(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two,first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five)))
            self.data[gene_name].update(pd.DataFrame({'element' : ['3ETS']}, index=range(first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five,first_flank+ETS+eighteen+its_one+five_dot_eight+its_two+twenty_five+ETS_end)))
            # self.data[gene_name].to_csv("olo_file.txt", sep='\t')

            data = self.data[gene_name]
            RDN37_data = data[data.element != 'unknown']

            # corr_dict = {"p" : "pearson" , "k" : "kendall" , "s" : "spearman"}
            corr_dict = {"p" : "pearson" , "s" : "spearman"}
            if output == 'a':
                print "# Calculating all correlations..."
                for i in corr_dict:
                    print "# Calculating correlations("+corr_dict[i]+")..."
                    matrix = RDN37_data.corr(method=corr_dict[i],min_periods=1)
                    matrix.to_csv(self.prefix+"RDN37_correlation_"+corr_dict[i]+".table", sep='\t')

            #calculate Pearson for different types
                    if elements == True:
                        for this_type in ['5ETS', '18S', 'ITS1', '5.8S', 'ITS2', '25S', '3ETS']:
                            new_data = RDN37_data[RDN37_data.element == this_type]
                            # #element internal ntotal normalization
                            # ntotal_data = pd.DataFrame()
                            # for exp_old in self.experiments:
                            #     exp = exp_old+'_ntotal'
                            #     ntotal_data[exp] = new_data[exp_old]/new_data[exp_old].sum()
                            # self.data[gene_name][exp] = self.data[gene_name][exp].add(0.000001) #adding pseudocounts

                            matrix = new_data.corr(method=corr_dict[i],min_periods=1)
                            matrix.to_csv(self.prefix+this_type+"_correlation_"+corr_dict[i]+".table", sep='\t')
                            # matrix_n = ntotal_data.corr(method=corr_dict[i],min_periods=1)
                            # matrix_n.to_csv(self.prefix+this_type+"_correlation_ntotal_"+corr_dict[i]+".table", sep='\t')
        return True

    def test_print(self, what):
        pd.set_option('display.max_rows', 5000)
        print what
        pd.reset_option('display.max_rows')
        return True