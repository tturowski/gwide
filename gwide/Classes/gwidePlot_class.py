#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 20110"
__version__		= "1.0"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

import numpy as np
import sys, collections, re
from pypeaks import Data
import pandas as pd
from pyCRAC.Parsers import GTF2
import matplotlib.pyplot as plt
import matplotlib
import plotly.plotly as py
import plotly.graph_objs as go
import scipy.stats.morestats as ss
import scipy.stats.stats as sss

class GenomeWidePlot():
    def __init__(self, gtf_file, five_prime_flank, three_prime_flank, hits_threshold, lookahead, prefix, readthrough_start, normalized, publish):
        self.gtf_file           =   str(gtf_file)
        self.gtf                =   GTF2.Parse_GTF()
        self.gtf.read_GTF(self.gtf_file)
        self.genes              =   dict()
        self.data               =   dict()
        self.id_to_names        =   dict()
        self.rt                 =   dict()                  # designated to work with one experiment only
        self.genes_name_list    =   list()
        self.genes_id_list      =   list()
        self.five_prime_flank   =   five_prime_flank
        self.three_prime_flank  =   three_prime_flank
        self.five_prime_to_print=   100
        self.longest_gene       =   0
        self.hits_threshold     =   hits_threshold
        self.lookahead          =   lookahead
        self.prefix             =   str(prefix)
        self.readthrough_start  =   readthrough_start
        self.experiments        =   list()
        self.normalized         =   normalized  # -n option allows for using normalized dataset (reads p[er Million)
        self.list_of_peaks = dict()
        self.publish            =   publish

    def read_csv(self, concat_file, skip_nucleotide=False):
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
            if self.longest_gene < gene_length:
                self.longest_gene = gene_length
            self.genes[gene_name] = {
                'gene_name'     :   gene_name,
                'gene_id'       :   gene_id,
                'gene_length'   :   gene_length,
                'introns'       :   self.get_introns(gene_name),   # [[len,len...][(start,stop),etc]] <-[[list_with_lengths],[list_with_start and stop]]
                'RT'            :   dict(),
                'total'         :   dict(),
                'total_av'      :   dict(),
                'total_med'     :   dict(),
                'total_std'     :   dict(),
                'a'             :   dict(),
                'b'             :   dict(),
                'i'             :   dict(),
                'e'             :   dict(),
                'f'             :   dict(),
                'd'             :   dict(),
                'd_av'          :   dict(),
                'd_med'         :   dict(),
                'd_std'         :   dict(),
                'g'             :   dict(),
                'g_av'          :   dict(),
                'g_med'         :   dict(),
                'g_std'         :   dict(),
                                    }
            self.genes[gene_name]['exons'] = self.get_exons(gene_name)
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
            if skip_nucleotide == False:
                self.data[gene_name]['nucleotide'] = gene_csv['nucleotide'][:(self.five_prime_flank + gene_length + self.three_prime_flank):]
            for e in self.experiments:
                self.data[gene_name][e] = gene_csv[gene_csv.exp_name == e]['hits']
                if self.normalized == True:
                    self.data[gene_name][(e+"_nrpm")] = gene_csv[gene_csv.exp_name == e]['n_hits']    #to work with data normalized reads per M
            self.data[gene_name] = self.data[gene_name].fillna(0)
        self.genes_id_list.sort()
        return True

    def get_introns(self, gene_name):
        gene_coord = self.gtf.chromosomeCoordinates(gene_name)
        introns_coord_raw = self.gtf.intronCoordinates(gene_name)
        if not introns_coord_raw:
            return [[],[]]
        else:
            introns = [[],[]]
            for i in range(0,len(introns_coord_raw)):
                intron_coord = list(introns_coord_raw[i])
                intron_len = max(intron_coord) - min(intron_coord) + 2 #corrected
                introns[0].append(intron_len)
                if self.gtf.strand(gene_name) == "+":
                    intron_start = min(intron_coord) - min(gene_coord)
                elif self.gtf.strand(gene_name) == "-":
                    intron_start = max(gene_coord) - max(intron_coord)
                intron_stop = intron_start + intron_len - 1 #corrected
                introns[1].append((intron_start,intron_stop))
            return introns

    def get_exons(self, gene_name):
        exons = list()
        introns = self.genes[gene_name]['introns']
        if len(introns[0]) == 0:
            exons = [[1, self.genes[gene_name]['gene_length']]]
            # print 'no introns for this gene'
        elif len(introns[0]) == 1:
            exons = [[1,introns[1][0][0]-1],[introns[1][0][1]+1, self.genes[gene_name]['gene_length']]]
        else:
            print "ERROR: Genes with more than one intron are not supported in this version."
        return exons

    def find_peaks(self):
        print '# Finding peaks...'
        for i in self.data:

            self.list_of_peaks[i] = dict()
            for e in self.experiments:
                self.list_of_peaks[i][e] = dict()
                hist = Data(list(self.data[i].index), list(self.data[i][e]), smoothness=1, default_smooth=False)
                hist.normalize()
                try:
                    hist.get_peaks(method="slope", peak_amp_thresh=0.000010, valley_thresh=0.00003, intervals=None,
                                   lookahead=self.lookahead, avg_interval=100)
                    self.list_of_peaks[i][e]['peaks'] = sorted(np.array(hist.peaks['peaks'][0]).tolist())
                    self.list_of_peaks[i][e]['valleys'] = sorted(np.array(hist.peaks['valleys'][0]).tolist())
                except ValueError:
                    pass
                    self.list_of_peaks[i][e]['peaks'] = []
                    self.list_of_peaks[i][e]['valleys'] = []
        # print self.list_of_peaks
        return True

    def calculate(self, details=True, ntotal=False, nmax=False, pscounts=False):
        if details == True or nmax == True or ntotal == True:
            print '# Calculating readthrough and others...'
        else:
            print '# Calculating readthrough...'

        for gene_name in self.genes:
            transcription_start =   self.five_prime_flank-21
            gene_length         =   self.genes[gene_name]['gene_length']
            gene_end            =   self.five_prime_flank + gene_length
            RT_begin            =   self.five_prime_flank + gene_length + self.readthrough_start
            gene_middle         =   self.five_prime_flank + ( gene_length / 2 )
            three_middle           =   self.five_prime_flank + gene_length + (self.three_prime_flank / 2)
            three_one_third        =   self.five_prime_flank + gene_length + (self.three_prime_flank / 3)
            three_two_third        =   self.five_prime_flank + gene_length + (2*(self.three_prime_flank / 3))

        #getting intron length
            introns = list()
            if not self.genes[gene_name]['introns'][0]:
                intron_length = 0
            else:
                for intron in range(0,len(self.genes[gene_name]['introns'][0])):
                    introns.append(str(self.genes[gene_name]['introns'][0][intron]))
                intron_length = int(''.join(map(str,introns)))
                intron_start_stop = self.genes[gene_name]['introns'][1][0]      #start and stop of first intron only!

            for exp in self.experiments:
        # !! changes in exp name !!
                exp_old = exp
                try:
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

                #slicing dataframes
                        total       =   self.data[gene_name][transcription_start:].sum()[exp]
                        total_av    =   self.data[gene_name][transcription_start:].mean()[exp]
                        total_med   =   self.data[gene_name][transcription_start:].median()[exp]
                        total_SD    =   self.data[gene_name][transcription_start:].std()[exp]
                        c           =   float(self.data[gene_name][RT_begin:].sum()[exp])
                        if details==True:
                            g       =   self.data[gene_name][transcription_start:gene_end].sum()[exp]
                            g_av    =   self.data[gene_name][transcription_start:gene_end].mean()[exp]
                            g_med   =   self.data[gene_name][transcription_start:gene_end].median()[exp]
                            g_SD    =   self.data[gene_name][transcription_start:gene_end].std()[exp]
                            a1      =   float(self.data[gene_name][transcription_start:gene_middle].sum()[exp])
                            a2      =   self.data[gene_name][gene_middle+1:gene_end].sum()[exp]
                            d       =   self.data[gene_name][gene_end+1:].sum()[exp]
                            d_av    =   self.data[gene_name][gene_end+1:].mean()[exp]
                            d_med   =   self.data[gene_name][gene_end+1:].median()[exp]
                            d_SD    =   self.data[gene_name][gene_end+1:].std()[exp]
                            e1      =   float(self.data[gene_name][gene_end+1:three_middle].sum()[exp])
                            e2      =   self.data[gene_name][three_middle+1:].sum()[exp]
                            f1      =   float(self.data[gene_name][gene_end+1:three_one_third].sum()[exp])
                            f2      =   self.data[gene_name][three_one_third+1:three_two_third].sum()[exp]
                            f3      =   self.data[gene_name][three_two_third+1:].sum()[exp]
                            if intron_length > 0:
                                b1  =   float(self.data[gene_name][transcription_start:intron_start_stop[0]+self.five_prime_flank].sum()[exp])
                                b2  =   float(self.data[gene_name][intron_start_stop[0]+1+self.five_prime_flank:intron_start_stop[1]+self.five_prime_flank].sum()[exp])
                                b3  =   self.data[gene_name][intron_start_stop[1]+1+self.five_prime_flank:gene_end].sum()[exp]
                #calculating
                        RT = np.float64(c) / total #allows for dividing by 0
                        if details==True:
                            a = np.float64(a1) / a2
                            e = np.float64(e1) / e2
                            f = np.float64(f1) / f3

                            if intron_length > 0:
                                b = np.float64(b1) / b3
                                i = np.float64(b2) / (b1 + b3)
                            else:
                                b = 0
                                i = 0
                    else:
                        RT, total, total_av, total_med, total_SD, c = ['too_low_reads'] * 6
                        if details == True:
                            g, g_av, g_med, g_SD, d, d_av, d_med, d_SD, a, e, f, b, i = ['too_low_reads'] * 13
                except KeyError as e:
                    # print "Error raised by key: "+str(e)
                    RT, total, total_av, total_med, total_SD, c = ['no_reads'] * 6
                    if details == True:
                            g, g_av, g_med, g_SD, d, d_av, d_med, d_SD, a, e, f, b, i = ['no_reads'] * 13
                #stroing in dictionary
                self.genes[gene_name]['RT'][exp_old] = RT
                if details==True:
                    self.genes[gene_name]['total'][exp_old] = total
                    self.genes[gene_name]['total_av'][exp_old]  = total_av
                    self.genes[gene_name]['total_med'][exp_old] = total_med
                    self.genes[gene_name]['total_std'][exp_old] = total_SD
                    self.genes[gene_name]['a'][exp_old]         = a
                    self.genes[gene_name]['b'][exp_old]         = b
                    self.genes[gene_name]['i'][exp_old]         = i
                    self.genes[gene_name]['e'][exp_old]         = e
                    self.genes[gene_name]['f'][exp_old]         = f
                    self.genes[gene_name]['d'][exp_old]     = d
                    self.genes[gene_name]['d_av'][exp_old]  = d_av
                    self.genes[gene_name]['d_med'][exp_old] = d_med
                    self.genes[gene_name]['d_std'][exp_old] = d_SD
                    self.genes[gene_name]['g'][exp_old]     = g
                    self.genes[gene_name]['g_av'][exp_old]  = g_av
                    self.genes[gene_name]['g_med'][exp_old] = g_med
                    self.genes[gene_name]['g_std'][exp_old] = g_SD

        ### below part of function is needed to sort tRNA according to readthrough
                tRNA_group = self.genes[gene_name]['gene_id'][0:2]
                gene_id = self.genes[gene_name]['gene_id']
                if tRNA_group not in self.rt:
                    self.rt[tRNA_group] = dict()
                self.rt[tRNA_group][gene_id] = RT
        return True

    def read_list(self, file):
        list_file = open(file,"r")
        self.to_include_list = list()
        self.align_position = dict()
        for line in list_file:
            if line[0] == "#": continue
            line_elements = line.strip().split('\t')
            gene_name, position = line_elements[0],int(line_elements[1])
            self.to_include_list.append(gene_name)
            self.align_position[gene_name] = position
        return self.to_include_list, self.align_position

    def filter_out(self, gene_name, filter, experiment_to_filter, current_exp):
        #function to filter out not matching genes
        if not filter:
            return True

        filter_elements = filter.strip().split('_')
        factor, compare, value = filter_elements[0], filter_elements[1], float(filter_elements[2])
        if factor not in ['RT', 'a', 'b', 'i', 'e', 'f', 'intron', 'introns']:
            exit("I can not recognize factor: "+str(factor))
        if compare not in ['a', 'b', 'above', 'below', '>', '<']:
            exit("I can not recognize comparator: "+str(factor))

        #for introns
        if factor == 'intron' or factor == 'introns' or factor == 'i':
            if len(self.genes[gene_name]['introns'][0]) > value:
                if compare in ['above', 'a', '>']:
                    return True
                else:
                    return False
            elif len(self.genes[gene_name]['introns'][0]) < value:
                if compare in ['below', 'b', '<']:
                    return True
                else:
                    return False
        else:
            #for rest of factors
            if experiment_to_filter:
                exp = experiment_to_filter
            else:
                exp = current_exp
            if self.genes[gene_name][factor][exp] > value:
                if compare in ['above', 'a', '>']:
                    return True
                else:
                    return False
            elif self.genes[gene_name][factor][exp] <= value:
                if compare in ['below', 'b', '<']:
                    return True
                else:
                    return False
            else:
                print "Mistake in -f statement: "+factor+compare+str(value)

    def filter_to_print(self, filter, experiment):
        #make filter to print as a text when plotting
        if not filter:
            return str()
        filter_elements = filter.strip().split('_')
        factor, compare, value = filter_elements[0], filter_elements[1], filter_elements[2]
        if compare in ['above', 'a', '>']:
            compare = 'above'
        elif compare in ['below', 'b', '<']:
            compare = 'below or equal'
        if experiment:
            return "Using filter: "+factor+" "+compare+" "+str(value)+" in dataset "+experiment
        else:
            return "Using filter: "+factor+" "+compare+" "+str(value)

    def categorize_tRNA(self,list_of_genes, output='text'):
        isoacceptors = list()
        anticodones = list()
        for gene_name in list_of_genes:
            gene_id = self.genes[gene_name]['gene_id']
            tRNA_isoacceptor = self.genes[gene_name]['gene_id'][0:2]
            tRNA_anticodon = self.genes[gene_name]['gene_id'][3:6]
            isoacceptors.append(tRNA_isoacceptor)
            anticodones.append(tRNA_anticodon)
        isoacceptors = list(set(isoacceptors))
        anticodones = list(set(anticodones))
        if output=='text':
            return 'representing '+str(len(isoacceptors))+' isotype(s) and '+str(len(anticodones))+' isoacceptor(s)'

    def group_experiments(self, a, b, exp_to_use):
        print 'Looking for groups of experiments...'
        ### function identifiying 'to_divide' and 'divisor' for ratio results plotting.
        ### Please use nomenclature gene_variant
        ### IMPORTANT: only two experiments may contain the same root(gene) added to 'a' or 'b'
        ### i.e. list_of_experiments = ['A190_MD', 'A190_total', 'A1310_MD', 'A1310_total'] additional A190 or A1310 is forbidden
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
            print "Unequal no. of experiments with -a and -b parameter:"
            print "-a exp: "+','.join(a_experiments)+"\n -b exp: "+','.join(b_experiments)
            exit()
        a_experiments.sort()
        b_experiments.sort()
        ratio_exp_list = list()
        if len(a_experiments) > 1 or len(b_experiments) > 1:
            for e in a_experiments:
                gene = re.search(r'(.*)_'+a, e).group(1)
                second = [d for d in b_experiments if re.search(gene,d)]
                if len(second) > 1:
                    print 'To many elements to divide for '+gene
                    exit()
                couple = [e+exp_to_use, second[0]+exp_to_use]
                ratio_exp_list.append(couple)
        else:
            ratio_exp_list.append([a_experiments[0]+exp_to_use,b_experiments[0]+exp_to_use])
        print 'Found pairs:'
        print ratio_exp_list
        return ratio_exp_list

    def plotSubplot(self, fig, layout, plot_no, title, data, line_color, select=None):
        fig.add_subplot(layout[0],layout[1],plot_no)
        fig.tight_layout()
        if not select:
            plt.title(title)
            plt.plot(data.index, data["sum"], color="black")
        elif select:
            plt.title(plot_no)
            data = data[select[0]:select[1]]
            data['zero'] = np.nan
            data['zero'] = data['zero'].fillna(0)
            plt.plot(data.index, data["sum"], color="black")
            plt.fill_between(data.index, data["sum"], data['zero'], color="grey", alpha=0.10)
        plt.axvline(self.five_prime_flank, color=line_color)
        plt.grid()
        return True

    def plotSinglePlot(self, title, data, line_color, online=False):
        print "Printing single plot in publication quality..."
        fig = plt.figure(figsize=(3, 2), facecolor='w', edgecolor='k')
        fig.add_subplot(1,1,1)
        plt.plot(data.index, data["sum"], color="black")
        plt.xlabel('position (nt)', fontsize=7)
        # plt.xticks(xxx)
        plt.ylabel('no. of reads', fontsize=7)
        plt.title(title, fontsize=9)
        plt.axvline(self.five_prime_flank, color=line_color)
        plt.tight_layout()
        plt.savefig(title+'.png', dpi=300)
        if online == True:
            py.iplot_mpl(fig, filename='plotly_')
        plt.clf()
        return True

    def plotlyFigure(self, trace, title, data_min, data_max, annotation_text=""):
        '''
        :param trace: pyplot trace object
        :param title: plot title
        :param annotation_text: text annotating line to align
        :return: pyplot figure object
        '''
        if annotation_text=="3' end":
            line_color="rgb(166, 28, 0)"
        elif annotation_text=="5' end":
            line_color="rgb(53, 118, 20)"
        else:
            line_color="rgb(0, 0, 0)"
        return {'data': [trace], 'layout': {'autosize':False,
                        'yaxis':{
                            'tickfont':{'color':"rgb(67, 67, 67)",'size':18},
                            'title':"no. of reads",
                            'range':[0,],
                           'titlefont':{'color':"rgb(67, 67, 67)"},
                            'type':"linear",
                            'autorange':True,
                            'exponentformat':"power"
                                },
                        'title':title,
                        'height':400,
                        'width':400,
                        'shapes': [{
                        # Line Vertical
                            'type': 'line',
                            'x0': 250,
                            'y0': data_min,
                            'x1': 250,
                            'y1': data_max,
                            'line': {'color': line_color,'width': 2}
                                }],
                        'xaxis':{
                            'tickfont':{'color':"rgb(67, 67, 67)",'size':18},
                            'title':"position (nt)",
                            'range':[0,],
                            'titlefont':{'color':"rgb(67, 67, 67)"},
                            'type':"linear",
                            'autorange':True,
                            'exponentformat':"power"
                                },
                        }}

    def plotlyDataToFigure(self, data, title, annotation_text=""):
        '''
        :param data: pandas dataframe object with data
        :param title: plot title
        :param annotation_text: text annotating line to align
        :return: pyplot figure object
        '''
        trace = go.Scatter(x=data.index, y=data['sum'], name = 'name', connectgaps=True)
        return self.plotlyFigure(trace=trace, title=title, annotation_text=annotation_text, data_min=min(data['sum']), data_max=max(data['sum']))

    def checkSelect(self,select):
        '''
        :param select: select parameter from option parser
        :return: ranges or syntax error
        '''
        ranges = select.strip().split('_')
        if len(ranges) > 2 or len(ranges) < 2:
            exit("--select - inapropriate no of elements, accepts only number1_number2")
        ranges = [int(i) for i in ranges]
        print "Plotting in ranges: "+str(ranges)
        return ranges

    def isEven(self, number):
        return number % 2 == 0

    def table(self, filter, experiment_to_filter):
        """saving *csv file with flanks 50 nt 5' and 250 nt 3' """
        for e in self.experiments:
            #initiating dataframes
            raw_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            no_of_genes = 0
            list_of_genes = list()
            for gene_name in self.genes_name_list:
                if self.filter_out(gene_name=gene_name, filter=filter, experiment_to_filter=experiment_to_filter, current_exp=e) == True:
                    no_of_genes += 1
                    list_of_genes.append(gene_name)
                # 5` aligned
                    raw_5data[gene_name] = self.data[gene_name][e]
            transposed = raw_5data.transpose()
            transposed = transposed.fillna(0)
            transposed.to_csv(str(e)+"_heatmap.csv")

    def ratio(self, to_divide, divisor, filter, exp_to_use=str(), select=None):
        new_exp_list = self.group_experiments(to_divide, divisor, exp_to_use=exp_to_use) #exp_to_use allows for normalizations
        if select:
            select = self.checkSelect(select)
        print "# Plotting genom wide plots..."
        #copying data
        for e in new_exp_list:
            #initiating dataframes (a for a exp, b for b exp and ratio)
            a_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            a_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            b_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            b_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            ratio_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            ratio_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            log2_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            log2_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            for gene_name in self.genes_name_list:
                gene_length = self.genes[gene_name]['gene_length']
                if self.filter_out(gene_name=gene_name, filter=filter, experiment_to_filter=e[1].strip('_ntotal'), current_exp=e[0]) == True:
                    a_5data[gene_name] = self.data[gene_name][e[0]]
                    temp_data_a = self.data[gene_name][e[0]][gene_length:].reset_index()
                    a_3data[gene_name] = temp_data_a[e[0]]

                    b_5data[gene_name] = self.data[gene_name][e[1]]
                    temp_data_b = self.data[gene_name][e[1]][gene_length:].reset_index()
                    b_3data[gene_name] = temp_data_b[e[1]]

                    ratio_5data[gene_name] = self.data[gene_name][e[0]]/self.data[gene_name][e[1]]
                    ratio_3data[gene_name] = temp_data_a[e[0]]/temp_data_b[e[1]]

                    log2_5data[gene_name] = np.log2(self.data[gene_name][e[0]]/self.data[gene_name][e[1]])
                    log2_3data[gene_name] = np.log2(temp_data_a[e[0]]/temp_data_b[e[1]])
            # sum
            a_5data['sum'] = a_5data.sum(axis=1)
            a_3data['sum'] = a_3data.sum(axis=1)
            b_5data['sum'] = b_5data.sum(axis=1)
            b_3data['sum'] = b_3data.sum(axis=1)
            ratio_5data['sum'] = ratio_5data.sum(axis=1)
            ratio_3data['sum'] = ratio_3data.sum(axis=1)
            log2_5data['sum'] = log2_5data.sum(axis=1)
            log2_3data['sum'] = log2_3data.sum(axis=1)

            # plotting
            fig = plt.figure(figsize=(14,9), facecolor='w', edgecolor='k')
            five = "5' aligned reads for "
            three = "3' aligned reads for "
            layout = [4,2]
            self.plotSubplot(fig=fig, layout=layout, plot_no=1, title=five+e[0], data=a_5data, line_color="green", select=select)
            self.plotSubplot(fig=fig, layout=layout, plot_no=2, title=three+e[0], data=a_3data, line_color="#7f0f0f", select=select)
            self.plotSubplot(fig=fig, layout=layout, plot_no=3, title=five+e[1], data=b_5data, line_color="green", select=select)
            self.plotSubplot(fig=fig, layout=layout, plot_no=4, title=three+e[1], data=b_3data, line_color="#7f0f0f", select=select)
            self.plotSubplot(fig=fig, layout=layout, plot_no=5, title=five+e[0]+"/"+e[1]+" ratio", data=ratio_5data, line_color="green", select=select)
            self.plotSubplot(fig=fig, layout=layout, plot_no=6, title=three+e[0]+"/"+e[1]+" ratio", data=ratio_3data, line_color="#7f0f0f", select=select)
            self.plotSubplot(fig=fig, layout=layout, plot_no=7, title=five+'log2 '+e[0]+"/"+e[1]+" ratio", data=log2_5data, line_color="green", select=select)
            self.plotSubplot(fig=fig, layout=layout, plot_no=8, title=three+'log2 '+e[0]+"/"+e[1]+" ratio", data=log2_3data, line_color="#7f0f0f", select=select)
            # fig.tight_layout()
            # font = {'family' : 'normal', 'weight' : 'normal', 'size'   : 20}
            # matplotlib.rc('font', **font)
            name = self.prefix+exp_to_use+to_divide+'_'+divisor
            if filter:
                name = name+'_'+filter
            if select:
                name = 'selected_'+name
            plt.savefig(name+'_ratio.png', dpi=300)
            plt.clf()
            
            if self.publish == True:
                a_5data.to_csv(self.prefix+exp_to_use+"a_5data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=a_5data, title=five+e[0], annotation_text="5' end"), filename=name+'_a5.png')
                a_3data.to_csv(self.prefix+exp_to_use+"a_3data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=a_3data, title=three+e[0], annotation_text="3' end"), filename=name+'_a3.png')
                b_5data.to_csv(self.prefix+exp_to_use+"b_5data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=b_5data, title=five+e[1], annotation_text="5' end"), filename=name+'_b5.png')
                b_3data.to_csv(self.prefix+exp_to_use+"b_3data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=b_3data, title=three+e[1], annotation_text="3' end"), filename=name+'_b3.png')
                ratio_5data.to_csv(self.prefix+exp_to_use+"ratio_5data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=ratio_5data, title=five+e[0]+"/"+e[1]+" ratio", annotation_text="5' end"), filename=name+'_ratio5.png')
                ratio_3data.to_csv(self.prefix+exp_to_use+"ratio_3data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=ratio_3data, title=three+e[0]+"/"+e[1]+" ratio", annotation_text="3' end"), filename=name+'_ratio3.png')
                log2_5data.to_csv(self.prefix+exp_to_use+"log2_5data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=log2_5data, title=five+'log2 '+e[0]+"/"+e[1]+" ratio", annotation_text="5' end"), filename=name+'_log2ratio5.png')
                log2_3data.to_csv(self.prefix+exp_to_use+"log2_3data.txt", sep="\t")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=log2_3data, title=three+'log2 '+e[0]+"/"+e[1]+" ratio", annotation_text="3' end"), filename=name+'_log2ratio3.png')
        return True

    def std(self, filter, experiment_to_filter, exp_to_use=str()):
        print "# Plotting genom wide plots..."
        #inintiating DataFrame for all experiments
        all_raw_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        all_raw_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        #copying data
        for e in self.experiments:
            #initiating dataframes
            raw_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            raw_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
            no_of_genes = 0
            list_of_genes = list()
            for gene_name in self.genes_name_list:
                if self.filter_out(gene_name=gene_name, filter=filter, experiment_to_filter=experiment_to_filter, current_exp=e+exp_to_use) == True:
                    no_of_genes += 1
                    list_of_genes.append(gene_name)
                # 5` aligned
                    raw_5data[gene_name] = self.data[gene_name][e+exp_to_use]
                # 3` aligned
                    gene_length =   self.genes[gene_name]['gene_length']
                    temp_data = self.data[gene_name][e+exp_to_use][gene_length:].reset_index()
                    raw_3data[gene_name] = temp_data[e+exp_to_use]

            # sum
            raw_5data['sum'] = raw_5data.sum(axis=1)
            raw_3data['sum'] = raw_3data.sum(axis=1)
            all_raw_5data[e+exp_to_use] = raw_5data.sum(axis=1)
            all_raw_3data[e+exp_to_use] = raw_3data.sum(axis=1)

            # text_file = open(self.prefix+"std"+e+".list", "w")
            # for i in list_of_genes:
            #     text_file.write(i + "\t" + self.genes[i]['gene_id'] + "\n")
            # text_file.close()

            # plotting 10` aligned
            print "# Printing plot for "+str(no_of_genes)+" genes "+self.categorize_tRNA(list_of_genes,output='text')+"..."
            matplotlib.use('Agg')
            fig = plt.figure(figsize=(12, 9), facecolor='w', edgecolor='k')
            five = '5` aligned raw reads for '
            three = '3` aligned raw reads for '
            layout = [3,2]
            self.plotSubplot(fig=fig, layout=layout, plot_no=1, title=five+e, data=raw_5data, line_color="green")
            self.plotSubplot(fig=fig, layout=layout, plot_no=2, title=three+e, data=raw_3data, line_color="#7f0f0f")
            fig.tight_layout()



            if filter == None:
                filter = str()
            name = self.prefix+exp_to_use+filter+e
            # save text file
            # text_file = open(name+".list", "w")
            # for i in list_of_genes:
            #     text_file.write(i + "\t" + self.genes[i]['gene_id'] + "\n")
            # text_file.close()
            #save .png file
            plt.savefig(name+'.png', dpi=200)
            plt.clf()

            if self.publish == True:
                raw_5data.to_csv(self.prefix+"raw_5data.txt", sep="\t")
                raw_3data.to_csv(self.prefix+"raw_3data.txt", sep="\t")
                # self.plotSinglePlot(title=five+e, data=raw_5data, line_color="green")
                # self.plotSinglePlot(title=three+e, data=raw_3data, line_color="#7f0f0f")
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=raw_5data, title=five+e, annotation_text="5' end"), filename=name+'_5.png')
                py.image.save_as(figure_or_data=self.plotlyDataToFigure(data=raw_3data, title=three+e, annotation_text="3' end"), filename=name+'_3.png')

                # #plot online
                # py.iplot(figure_or_data=self.plotlyDataToFigure(data=raw_5data, title=five+e, annotation_text="5' end"), filename=name+'_5.png')
                # py.iplot(figure_or_data=self.plotlyDataToFigure(data=raw_3data, title=three+e, annotation_text="3' end"), filename=name+'_3.png')
        all_raw_5data.to_csv(self.prefix+exp_to_use+"ALL_raw_5data.txt", sep="\t")
        all_raw_3data.to_csv(self.prefix+exp_to_use+"ALL_raw_3data.txt", sep="\t")
        return True

    def aligner(self, file, filter, experiment_to_filter):
        print "# Plotting genom wide plots with chosen aligner..."
        #initiating dataframes
        raw_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        raw_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        raw_chosen_data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])

        #copying data
        for e in self.experiments:
            no_of_genes = 0
            list_of_genes = list()
            for gene_name in self.to_include_list:
                if gene_name in self.genes_name_list:
                    if self.filter_out(gene_name=gene_name, filter=filter, experiment_to_filter=experiment_to_filter, current_exp=e) == True:
                        no_of_genes += 1
                        list_of_genes.append(gene_name)
                    # 10` aligned
                        raw_5data[gene_name] = self.data[gene_name][e]
                    # 3` aligned
                        gene_length =   self.genes[gene_name]['gene_length']
                        temp_data = self.data[gene_name][e][gene_length:].reset_index()
                        raw_3data[gene_name] = temp_data[e]
                    # chosen aligner
                        aligner = self.align_position[gene_name]
                        temp_data = self.data[gene_name][e][aligner:].reset_index()
                        raw_chosen_data[gene_name] = temp_data[e]

            # sum
            raw_5data['sum'] = raw_5data.sum(axis=1)
            raw_3data['sum'] = raw_3data.sum(axis=1)
            raw_chosen_data['sum'] = raw_chosen_data.sum(axis=1)

            # plotting
            print "# Printing plot for "+str(no_of_genes)+" genes "+self.categorize_tRNA(list_of_genes,output='text')+"..."
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            # plt.title(str(no_of_genes)+" tRNA genes with in "+e+" dataset. "+self.filter_to_print(filter=filter, experiment=experiment_to_filter)+"\n"
            #           +self.categorize_tRNA(list_of_genes,output='text'),y=1.04)

            fig.add_subplot(3,1,1)
            plt.title('10` aligned raw reads')
            plt.plot(raw_5data.index, raw_5data['sum'], color="black")
            plt.axvline(self.five_prime_flank, color="green")
            plt.grid()

            fig.add_subplot(3,1,2)
            plt.title('3` aligned raw reads')
            plt.plot(raw_3data.index, raw_3data['sum'], color="black")
            plt.axvline(self.five_prime_flank, color="#7f0f0f")
            plt.grid()

            fig.add_subplot(3,1,3)
            plt.title('raw reads aligned according to '+file)
            plt.plot(raw_chosen_data.index, raw_chosen_data['sum'], color="black")
            plt.axvline(self.five_prime_flank, color="black")
            plt.grid()

            if filter == None:
                filter = str()
            name = self.prefix+filter+e
            #save text file
            text_file = open(name+".list", "w")
            for i in list_of_genes:
                text_file.write(i + "\t" + self.genes[i]['gene_id'] + "\t" + str(self.align_position[i]) + "\n")
            text_file.close()
            #save .png file
            plt.savefig(name+'.png')
            plt.clf()
        return True

    def makeRTGTF(self):
        print "# Making GTF with RT coordinates"
        name = self.prefix+"RTextension.gtf"
        new_GTF_file = open(name, "w")
        for gene_name in self.genes_name_list:
            highest_RT = 0
            for e in self.experiments:
                if not self.list_of_peaks[gene_name][e]['valleys']: continue
                RT_length = max(self.list_of_peaks[gene_name][e]['valleys'])-self.five_prime_flank
                # print "RT for "+gene_name+" in "+e+": "+str(RT_length)
                if RT_length > highest_RT: highest_RT = RT_length
            # print gene_name +' '+ str(highest_RT)
            # print self.gtf.chromosomeCoordinates(gene_name)
            # print self.gtf.genes[gene_name]['strand']
            if self.gtf.genes[gene_name]['strand'] == '+': start, stop = self.gtf.chromosomeCoordinates(gene_name)[1]+1, self.gtf.chromosomeCoordinates(gene_name)[0]+1+highest_RT
            elif self.gtf.genes[gene_name]['strand'] == '-': start, stop = self.gtf.chromosomeCoordinates(gene_name)[1]-highest_RT, self.gtf.chromosomeCoordinates(gene_name)[0]
            line = self.gtf.genes[gene_name]['chromosome']+'\t'\
                   +'tRNAextension'+"\t"\
                   +'exon'+'\t'\
                   +str(start)+'\t'\
                   +str(stop)+"\t.\t"\
                   +self.gtf.genes[gene_name]['strand']+"\t.\t"\
                   +'gene_id "'+self.gtf.genes[gene_name]['gene_id']+'"; transcript_id "'+self.gtf.genes[gene_name]['gene_id']\
                   +'"; gene_name "'+gene_name+'"; gene_biotype "tRNAextension"; transcript_name "'+gene_name+'";'
            if highest_RT != 0: new_GTF_file.write(line+'\n')
        new_GTF_file.close()

    def maketranscriptGTF(self):
        print "# Making GTF with transcript coordinates"
        print "# TSS is setup to -11"
        name = self.prefix+"transcript.gtf"
        new_GTF_file = open(name, "w")
        for gene_name in self.genes_name_list:
            highest_RT = 0
            for e in self.experiments:
                if not self.list_of_peaks[gene_name][e]['valleys']:
                    highest_RT = self.genes[gene_name]['gene_length']
                else:
                    RT_length = max(self.list_of_peaks[gene_name][e]['valleys'])-self.five_prime_flank
                    if RT_length > highest_RT: highest_RT = RT_length
            if self.gtf.genes[gene_name]['strand'] == '+': start, stop = self.gtf.chromosomeCoordinates(gene_name)[0]-10, self.gtf.chromosomeCoordinates(gene_name)[0]+1+highest_RT
            elif self.gtf.genes[gene_name]['strand'] == '-': start, stop = self.gtf.chromosomeCoordinates(gene_name)[1]-highest_RT, self.gtf.chromosomeCoordinates(gene_name)[1]+11
            line = self.gtf.genes[gene_name]['chromosome']+'\t'\
                   +'tRNA'+"\t"\
                   +'exon'+'\t'\
                   +str(start)+'\t'\
                   +str(stop)+"\t.\t"\
                   +self.gtf.genes[gene_name]['strand']+"\t.\t"\
                   +'gene_id "'+self.gtf.genes[gene_name]['gene_id']+'"; transcript_id "'+self.gtf.genes[gene_name]['gene_id']\
                   +'"; gene_name "'+gene_name+'"; gene_biotype "tRNA"; transcript_name "'+gene_name+'";'
            new_GTF_file.write(line+'\n')
        new_GTF_file.close()

    def printTrancriptLength(self):
        print "# Making GTF with transcript coordinates"
        print "# TSS is setup to -11"
        name = self.prefix+"transcript_len.txt"
        data = pd.DataFrame()
        for gene_name in self.genes_name_list:
            for e in self.experiments:
                if not self.list_of_peaks[gene_name][e]['valleys']:
                    transcript_len = self.genes[gene_name]['gene_length']
                else:
                    transcript_len = self.genes[gene_name]['gene_length'] + max(self.list_of_peaks[gene_name][e]['valleys']) - self.five_prime_flank + 11
                data.loc[gene_name, 'gene_id'] = self.gtf.genes[gene_name]['gene_id']
                data.loc[gene_name, e] = transcript_len
        data.to_csv(name, sep='\t')

    def Tdensity(self, peak_min, size):
        if not self.isEven(size):
            exit("Size should be even number.")
        print "# Calculating genome-wide T-density among peaks and valleys"
        peaks_T = list()
        valleys_T = list()
        for e in self.experiments:
            for gene_name in self.genes_name_list:
                # if self.filter_out(gene_name=gene_name, filter='RT_a_0.0001', experiment_to_filter=e, current_exp=e) == True:
                if not self.list_of_peaks[gene_name][e]['valleys'] or not self.list_of_peaks[gene_name][e]['peaks']: continue
                self.list_of_peaks[gene_name][e]['RT_peaks'] = list()
                self.list_of_peaks[gene_name][e]['RT_valleys'] = list()
                gene_length = self.genes[gene_name]['gene_length']
                for peak in self.list_of_peaks[gene_name][e]['peaks']:
                    if peak > gene_length+self.five_prime_flank:
                        next_valley = min([v for v in self.list_of_peaks[gene_name][e]['valleys'] if v > peak])
                        peak_T_content = self.data[gene_name][peak-(size/2):peak+(size/2)][self.data[gene_name].nucleotide == 'T']
                        next_valley_T_content = self.data[gene_name][next_valley-(size/2):next_valley+(size/2)][self.data[gene_name].nucleotide == 'T']
                        if peak_T_content.mean()[e] > peak_min:
                            peaks_T.append(float(peak_T_content.count()['nucleotide'])/size)
                            valleys_T.append(float(next_valley_T_content.count()['nucleotide'])/size)
        print 'Minimum peak average setup to '+str(peak_min)
        print 'Calculating for '+str(size)+" nt surrounding each peak and through(valley)"
        print 'No. of peaks: '+str(len(peaks_T))
        print 'No. of valleys: '+str(len(valleys_T))
        df = pd.DataFrame({'peaks_T' : peaks_T, 'valleys_T' : valleys_T})
        print 'Mean:'
        print df.mean()
        # print 'Madian:'
        # print df.median()
        # print 'SD:'
        # print df.std()
        print 'Wilcoxon:'
        wilcoxon = ss.wilcoxon(peaks_T, valleys_T)
        print wilcoxon
        print 'T-test (related samples):'
        # print sss.ttest_ind(peaks_T, valleys_T, equal_var=True)
        # print sss.ttest_ind(peaks_T, valleys_T, equal_var=False)
        print sss.ttest_rel(peaks_T, valleys_T)
        print 'To_copy: '+str(df.mean()['peaks_T'])+'\t'+str(df.mean()['valleys_T'])+'\t'+str(wilcoxon[1])

        df.to_csv('T_in_peaks_valleys.txt', sep='\t')
        # # Create a figure instance
        # fig = plt.figure(1, figsize=(9, 6))
        #         # Create an axes instance
        # ax = fig.add_subplot(111)
        #         # Create the boxplot
        # ax.set_xticklabels(['peaks', 'valleys'])
        # ax.get_xaxis().tick_bottom()
        # bp = ax.boxplot([peaks_T, valleys_T], notch=True, sym=None, vert=True, whis=[110,810],
        #         positions=None, widths=None, patch_artist=True,
        #         bootstrap=None, usermedians=None, conf_intervals=None,
        #         meanline=True, showmeans=False, showcaps=True,
        #         showbox=True, showfliers=True, boxprops=None, labels=None,
        #         flierprops=None, medianprops=None, meanprops=None,
        #         capprops=None, whiskerprops=None, manage_xticks=True)
        #
        # ## change outline color, fill color and linewidth of the boxes
        # for box in bp['boxes']:
        #     # change outline color
        #     box.set( color='black', linewidth=2)
        #     # change fill color
        #     box.set( facecolor ='grey')
        #
        # ## change color and linewidth of the whiskers
        # for whisker in bp['whiskers']:
        #     whisker.set(color='grey', linewidth=2)
        #
        # ## change color and linewidth of the caps
        # for cap in bp['caps']:
        #     cap.set(color='grey', linewidth=2)
        #
        # ## change color and linewidth of the medians
        # for median in bp['medians']:
        #     median.set(color='black', linewidth=1)
        #
        # ## change the style of fliers and their fill
        # for flier in bp['fliers']:
        #     flier.set(marker='o', color='grey', alpha=0.10)
        #
        # # Save the figure
        # name = self.prefix+'boxplot.png'
        # fig.savefig(name, bbox_inches='tight')

    def RT_aligner(self, filter, experiment_to_align):
        print "# Plotting genom wide plots with chosen aligner..."
        #initiating dataframes
        raw_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        raw_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        raw_RT3_data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])

        #copying data
        for e in self.experiments:
            if filter == None:
                filter = str()
            name = self.prefix+'RTaligned_'+filter+e

            if experiment_to_align:
                e_for_RTend = experiment_to_align
            else:
                e_for_RTend = e

            no_of_genes = 0
            list_of_genes = list()
            for gene_name in self.genes_name_list:
                if self.filter_out(gene_name=gene_name, filter=filter, experiment_to_filter=experiment_to_align, current_exp=e) == True:
                    no_of_genes += 1
                    list_of_genes.append(gene_name)
                # 10` aligned
                    raw_5data[gene_name] = self.data[gene_name][e]
                # 3` aligned
                    gene_length =   self.genes[gene_name]['gene_length']
                    temp_data = self.data[gene_name][e][gene_length:].reset_index()
                    raw_3data[gene_name] = temp_data[e]
                # chosen aligner
                    if not self.list_of_peaks[gene_name][e_for_RTend]['valleys']: continue
                    RT_length = max(self.list_of_peaks[gene_name][e_for_RTend]['valleys'])-self.five_prime_flank
                    temp_data = self.data[gene_name][e][RT_length:].reset_index()
                    raw_RT3_data[gene_name] = temp_data[e]

           #save text file
            text_file = open(name+".list", "w")
            for i in list_of_genes:
                if not self.list_of_peaks[i][e_for_RTend]['valleys']:
                    text_file.write(i + "\t" + self.genes[i]['gene_id'] + "\t0"+ "\n")
                else:
                    text_file.write(i + "\t" + self.genes[i]['gene_id'] + "\t" + str(max(self.list_of_peaks[i][e_for_RTend]['valleys'])-self.five_prime_flank) + "\n")
            text_file.close()

            # sum
            raw_5data['sum'] = raw_5data.sum(axis=1)
            raw_3data['sum'] = raw_3data.sum(axis=1)
            raw_RT3_data['sum'] = raw_RT3_data.sum(axis=1)

            # plotting
            print "# Printing plot for "+str(no_of_genes)+" genes "+self.categorize_tRNA(list_of_genes,output='text')+"..."
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            plt.title(str(no_of_genes)+" tRNA genes with in "+e+" dataset. "+self.filter_to_print(filter=filter, experiment=experiment_to_align)+"\n"
                      +self.categorize_tRNA(list_of_genes,output='text'),y=1.04)

            fig.add_subplot(3,1,1)
            plt.title('10` aligned raw reads')
            plt.plot(raw_5data.index, raw_5data['sum'])
            plt.axvline(self.five_prime_flank, color="green")

            fig.add_subplot(3,1,2)
            plt.title('3` aligned raw reads')
            plt.plot(raw_3data.index, raw_3data['sum'])
            plt.axvline(self.five_prime_flank, color="green")

            fig.add_subplot(3,1,3)
            if experiment_to_align:
                plt.title("raw reads aligned according to 3` end of read-through in exp: "+experiment_to_align)
            else:
                plt.title("raw reads aligned according to 3` end of read-through")
            plt.plot(raw_RT3_data.index, raw_RT3_data['sum'])
            plt.axvline(self.five_prime_flank, color="green")

            #save .png file
            plt.savefig(name+'.png')
            plt.clf()
        return True

    def genome_wide_plots_3end_normalizations(self, set_up=None, factor=None, lower_than=False):
        """examples of use:
        analyse_tRNA_from_concatv2.py -c glu_wt_C2.concat -o gwideToolkit -f RT -s 0.110
        analyse_tRNA_from_concatv2.py -c glu_wt_C2.concat -o gwideToolkit -f RT -s 0.110 --lower_than"""

        print "# Plotting genom wide plots..."
        #initiating dataframes
        raw_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        nmax_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        ntotal_5data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        raw_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        nmax_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        ntotal_3data = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        raw_3RTdata = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        nmax_3RTdata = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])
        ntotal_3RTdata = pd.DataFrame(index=range(1,(self.five_prime_flank + self.longest_gene + self.three_prime_flank+1)), columns=[])

        #copying data
        for e in self.experiments:
            no_of_genes = 0
            list_of_genes = list()
            for gene_name in self.genes_name_list:
                if self.list_of_peaks[gene_name][e]['valleys']:
                    no_of_genes += 1
                    list_of_genes.append(gene_name)
                # 10` aligned
                    raw_5data[gene_name] = self.data[gene_name][e]
                    nmax_5data[gene_name] = self.data[gene_name][e+'_nmax']
                    ntotal_5data[gene_name] = self.data[gene_name][e+'_ntotal']
                # 3` aligned
                    gene_length =   self.genes[gene_name]['gene_length']
                    temp_data = self.data[gene_name][e][gene_length:].reset_index()
                    raw_3data[gene_name] = temp_data[e]
                    temp_data = self.data[gene_name][e+'_nmax'][gene_length:].reset_index()
                    nmax_3data[gene_name] = temp_data[e+'_nmax']
                    temp_data = self.data[gene_name][e+'_ntotal'][gene_length:].reset_index()
                    ntotal_3data[gene_name] = temp_data[e+'_ntotal']
                # 3` od RT aligned
                    RT_length = max(self.list_of_peaks[gene_name][e]['valleys'])-self.five_prime_flank
                    temp_data = self.data[gene_name][e][RT_length:].reset_index()
                    raw_3RTdata[gene_name] = temp_data[e]
                    temp_data = self.data[gene_name][e+'_nmax'][RT_length:].reset_index()
                    nmax_3RTdata[gene_name] = temp_data[e+'_nmax']
                    temp_data = self.data[gene_name][e+'_ntotal'][RT_length:].reset_index()
                    ntotal_3RTdata[gene_name] = temp_data[e+'_ntotal']

            # sum
            raw_5data['sum'] = raw_5data.sum(axis=1)
            nmax_5data['sum'] = nmax_5data.sum(axis=1)
            ntotal_5data['sum'] = ntotal_5data.sum(axis=1)
            raw_3data['sum'] = raw_3data.sum(axis=1)
            nmax_3data['sum'] = nmax_3data.sum(axis=1)
            ntotal_3data['sum'] = ntotal_3data.sum(axis=1)
            raw_3RTdata['sum'] = raw_3RTdata.sum(axis=1)
            nmax_3RTdata['sum'] = nmax_3RTdata.sum(axis=1)
            ntotal_3RTdata['sum'] = ntotal_3RTdata.sum(axis=1)

            text_file = open(self.prefix+e+"_RTend.list", "w")
            for i in list_of_genes:
                text_file.write(i + "\t" + self.genes[i]['gene_id'] + "t" + str(max(self.list_of_peaks[i][e]['valleys'])-self.five_prime_flank) + "\n")
            text_file.close()

            # plotting 10` aligned
            print "# Printing plot for "+str(no_of_genes)+" genes..."
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            plt.title(str(no_of_genes)+" tRNA genes with in "+e+" dataset",y=1.04)
            fig.add_subplot(3,3,1)
            plt.title('10` aligned raw reads')
            plt.plot(raw_5data.index, raw_5data['sum'])
            plt.axvline(2100, color="green")
            fig.add_subplot(3,3,2)
            plt.title('3` aligned raw reads')
            plt.plot(raw_3data.index, raw_3data['sum'])
            plt.axvline(2100, color="green")
            fig.add_subplot(3,3,3)
            plt.title('3` of RT aligned raw reads')
            plt.plot(raw_3data.index, raw_3RTdata['sum'])
            plt.axvline(2100, color="green")
            fig.add_subplot(3,3,4)
            plt.title('10` aligned normalized (nmax)')
            plt.plot(nmax_5data.index, nmax_5data['sum'])
            fig.add_subplot(3,3,10)
            plt.title('3` aligned normalized (nmax)')
            plt.plot(nmax_3data.index, nmax_3data['sum'])
            fig.add_subplot(3,3,6)
            plt.title('3` of RT aligned normalized (nmax)')
            plt.plot(nmax_3data.index, nmax_3RTdata['sum'])
            fig.add_subplot(3,3,7)
            plt.title('10` aligned normalized (ntotal)')
            plt.plot(ntotal_5data.index, ntotal_5data['sum'])
            fig.add_subplot(3,3,8)
            plt.title('3` aligned normalized (ntotal)')
            plt.plot(ntotal_3data.index, ntotal_3data['sum'])
            fig.add_subplot(3,3,9)
            plt.title('3` of RT aligned normalized (ntotal)')
            plt.plot(ntotal_3data.index, ntotal_3RTdata['sum'])
            # fig.tight_layout()
            plt.savefig(self.prefix+e+'_RTaligned.png')
            plt.clf()
        return True