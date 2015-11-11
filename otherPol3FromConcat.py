#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "1.3"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

import numpy as np
import sys, collections, re
from pypeaks import Data
from pyCRAC.Parsers import GTF2
import matplotlib.pyplot as plt
import pandas as pd

#this class is slightly modified tRNAFromConcat

class OtherPol3FromConcat():
    def __init__(self, gtf_file, five_prime_flank, three_prime_flank, hits_threshold, lookahead, prefix, print_valleys, print_peaks, readthrough_start, normalized):
        self.gtf_file           =   str(gtf_file)
#        self.fasta_file         =   str(fasta_file)
        self.gtf                =   GTF2.Parse_GTF()
        self.gtf.read_GTF(self.gtf_file)
#        self.gtf.read_FASTA(self.fasta_file)
        self.genes              =   dict()
        self.data               =   dict()
        self.id_to_names        =   dict()
        self.rt                 =   dict()                  # designated to work with one experiment only
        self.genes_name_list    =   list()
        self.genes_id_list      =   list()
        self.five_prime_flank   =   five_prime_flank
        self.three_prime_flank  =   three_prime_flank
        self.hits_threshold     =   hits_threshold
        self.lookahead          =   lookahead
        self.prefix             =   str(prefix)
        self.print_peaks        =   print_peaks
        self.print_valleys      =   print_valleys
        self.readthrough_start  =   readthrough_start
        self.experiments        =   list()
        self.normalized         =   normalized  # -n option allows for using normalized dataset (reads p[er Million)
        if self.print_peaks == True or self.print_valleys == True:
            self.list_of_peaks = dict()

    def read_concat_file(self, concat_file, null_substitution):
        print '# Reading concat file.'
        for line in concat_file:
            try:
                if line[0] == "#": continue
                line_elements = line.strip().split('\t')
                if self.normalized == False:
                    gene_name, position, nucleotide, hits, experiment = line_elements[0],int(line_elements[1]), line_elements[2], int(line_elements[3]), line_elements[6]
                elif self.normalized == True:
                    gene_name, position, nucleotide, hits, experiment = line_elements[0],int(line_elements[1]), line_elements[2], int(line_elements[7]), line_elements[6]

        #adding new entry in genes dict
                if gene_name not in self.genes:
                    gene_id = self.gtf.genes[gene_name]['gene_id']
                    gene_length = self.gtf.geneLength(gene_name)
                    self.genes[gene_name] = {
                        'gene_name'     :   gene_name,
                        'gene_id'       :   gene_id,
                        'gene_length'   :   gene_length,
                        'introns'       :   self.get_introns(gene_name),   # [[len,len...][(start,stop),etc]] <-[[list_with_lengths],[list_with_start and stop]]
                        'RT'            :   dict()
                                            }
                    self.genes[gene_name]['exons'] = self.get_exons(gene_name)
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
                if position <= self.five_prime_flank:
                    self.data[gene_name][experiment].loc[position, 'position'] = position - self.five_prime_flank
                else:
                    self.data[gene_name][experiment].loc[position, 'position'] = position - self.five_prime_flank
                self.data[gene_name][experiment].loc[position, 'nucleotides'] = nucleotide
                if null_substitution == True and hits == 0: ## option necessry for ratio and log calculation
                    self.data[gene_name][experiment].loc[position, 'hits'] = 1
                else:
                    self.data[gene_name][experiment].loc[position, 'hits'] = hits

            except IndexError:
                sys.stderr.write("\nIndexError at line:\n%s\n" % line)
                pass
        self.genes_name_list.sort()
        self.genes_id_list.sort()
        self.experiments.sort()

    def get_introns(self, gene_name):
        gene_coord = self.gtf.chromosomeCoordinates(gene_name)
        introns_coord_raw = self.gtf.intronCoordinates(gene_name)
        if not introns_coord_raw:
            return [[],[]]
        else:
            introns = [[],[]]
            for i in range(0,len(introns_coord_raw)):
                intron_coord = list(introns_coord_raw[i])
                intron_len = max(intron_coord) - min(intron_coord)
                introns[0].append(intron_len)
                if self.gtf.strand(gene_name) == "+":
                    intron_start = min(intron_coord) - min(gene_coord)
                elif self.gtf.strand(gene_name) == "-":
                    intron_start = max(gene_coord) - max(intron_coord)
                intron_stop = intron_start + intron_len
                introns[1].append((intron_start,intron_stop))
            return introns

    def get_exons(self, gene_name):
        exons = list()
        introns = self.genes[gene_name]['introns']
        if len(introns[0]) == 0:
            exons = [[0, self.genes[gene_name]['gene_length']]]
            # print 'no introns for this gene'
        elif len(introns[0]) == 1:
            exons = [[0,introns[1][0][0]-1],[introns[1][0][1]+1, self.genes[gene_name]['gene_length']]]
        else:
            print "Genes with more than one intron are not supported in this version."
        return exons

    def find_peaks(self):
        print '# Finding peaks...'
        for i in self.data:
            self.list_of_peaks[i] = dict()
            for e in self.data[i]:
                self.list_of_peaks[i][e] = dict()
                hist = Data(list(self.data[i][e].index), list(self.data[i][e]['hits']), smoothness=1, default_smooth=False)
                hist.normalize()
                try:
                    hist.get_peaks(method="slope", peak_amp_thresh=0.00005, valley_thresh=0.00003, intervals=None,
                                   lookahead=self.lookahead, avg_interval=100)
                    self.list_of_peaks[i][e]['peaks'] = sorted(np.array(hist.peaks['peaks'][0]).tolist())
                    self.list_of_peaks[i][e]['valleys'] = sorted(np.array(hist.peaks['valleys'][0]).tolist())
                except ValueError:
                    pass
                    self.list_of_peaks[i][e]['peaks'] = []
                    self.list_of_peaks[i][e]['valleys'] = []
        return True

    def calculate_read_through(self):
        print '# Calculating readthrough.'
        for i in self.genes:
            RT_begin = self.five_prime_flank + self.genes[i]['gene_length'] + self.readthrough_start
            for e in self.experiments:
                gene_sum    =   float()
                RT_sum      =   float()
                try:
                    histogram   =   list(self.data[i][e]['hits'])
                    if max(histogram) >= self.hits_threshold:
                        for n in range((self.five_prime_flank-20),len(histogram)):
                            gene_sum += histogram[n]
                            if n > RT_begin:
                                RT_sum += histogram[n]
                            RT = RT_sum / gene_sum
                    else:
                        RT = 'too_low_reads'
                    self.genes[i]['RT'][e] = RT
                except KeyError:
                    RT = 'no_reads'
                    self.genes[i]['RT'][e] = RT
        ### below part of function is needed to sort tRNA according to readthrough
                tRNA_group = self.genes[i]['gene_id'][0:2]
                gene_id = self.genes[i]['gene_id']
                if tRNA_group not in self.rt:
                    self.rt[tRNA_group] = dict()
                self.rt[tRNA_group][gene_id] = RT
        return True

# making output; if -p option then text file, else to standard output
    def make_text_file(self, filename):
        print '# Making text file.'
        output_dict = dict()

        for i in self.genes:
            output_dict[i] = list()
            introns = list()
            output_dict[i].append(str(self.genes[i]['gene_name']))
            output_dict[i].append(str(self.genes[i]['gene_id']))
            output_dict[i].append(str(self.genes[i]['gene_length']))
            if not self.genes[i]['introns'][0]:
                introns_to_print = 'none'
            else:
                for intron in range(0,len(self.genes[i]['introns'][0])):
                    introns.append(str(self.genes[i]['introns'][0][intron]))
                    introns_to_print = ', '.join(introns)
            output_dict[i].append(introns_to_print) # introns length
            for e in self.experiments:
                output_dict[i].append(str(self.genes[i]['RT'][e]))
        experiments = '\t'.join(self.experiments)
        filename.write("# analyse_tRNA_pileups output file using gtf file:"+"\n")
        filename.write("# "+self.gtf_file+"\n")
        filename.write("# 5` flank: "+str(self.five_prime_flank)+"\n")
        filename.write("# 3` flank: "+str(self.three_prime_flank)+"\n")
        if self.normalized == True:
            filename.write("# calculations performed on data normalized to reads per Million \n")
        else:
            filename.write("# calculations performed on non-normalized data \n")
        filename.write("# readthrough calculated: last nucleotide of gene + "+str(self.readthrough_start)+" nt\n")
        filename.write("# threshold of hits: "+str(self.hits_threshold)+"\n")
        filename.write("# lookahead option for pypeaks: "+str(self.lookahead)+"\n")
        filename.write("# name"+"\t"+"ID"+"\t"+"length"+"\t"+"int_len"+"\t"+experiments+"\n")
        for g in self.genes_id_list:
            line = self.id_to_names[g]
            filename.write('\t'.join(output_dict[line]) + '\n')
        if filename != sys.stdout:
            filename.close()
        return True

### preparing dataframes to for plottingyeast_all_tRNA.list
    def create_features_to_plot(self):
        for i in self.data:
            frame_index = range(1,(self.five_prime_flank + self.genes[i]['gene_length'] + self.three_prime_flank+1))
            RT_begin = self.five_prime_flank + self.genes[i]['gene_length'] + self.readthrough_start
            for e in self.data[i]:
                histogram = list(self.data[i][e]['hits'])
                hist_max = max(histogram)
    # histogram for peaks and valleys
                if self.print_peaks == True or self.print_valleys == True:
                    peaks = self.list_of_peaks[i][e]['peaks']
                    valleys = self.list_of_peaks[i][e]['valleys']
                    for c in frame_index:
                        if c in peaks:
                            self.data[i][e].loc[c, 'peaks_to_plot'] = hist_max
                        else:
                            self.data[i][e].loc[c, 'peaks_to_plot'] = 0

                        if c in valleys:
                            self.data[i][e].loc[c, 'valleys_to_plot'] = hist_max
                        else:
                            self.data[i][e].loc[c, 'valleys_to_plot'] = 0

    # histogram for readthrough
                for c in frame_index:
                    if c > RT_begin:
                        self.data[i][e].loc[c, 'RT_to_plot'] = self.data[i][e].loc[c, 'hits']
                    else:
                        self.data[i][e].loc[c, 'RT_to_plot'] = 0
        return True

    def termination_efficency_valleys(self):
        print '# Calculating energy for termination efficiency...'
        dna_dna = {	'AA' : [7.9, 0.0222],
                    'TT' : [7.9, 0.0222],
                    'AT' : [7.2, 0.0204],
                    'TA' : [7.2, 0.0213],
                    'CA' : [8.5, 0.0227],
                    'TG' : [8.5, 0.0227],
                    'GT' : [8.4, 0.0224],
                    'AC' : [8.4, 0.0224],
                    'CT' : [7.8, 0.021],
                    'AG' : [7.8, 0.021],
                    'GA' : [8.2, 0.0222],
                    'TC' : [8.2, 0.0222],
                    'CG' : [10.6, 0.0272],
                    'GC' : [9.8, 0.0244],
                    'GG' : [8.0, 0.0199],
                    'CC' : [8.0, 0.0199]
                        }
        rna_dna = { 'AA' : [7.8, 0.0219],
                    'TT' : [11.5, 0.0364],
                    'AT' : [8.3, 0.0239],
                    'TA' : [7.8, 0.0232],
                    'CA' : [9.0, 0.0261],
                    'TG' : [10.4, 0.0284],
                    'GT' : [7.8, 0.0216],
                    'AC' : [5.9, 0.0123],
                    'CT' : [7.0, 0.0197],
                    'AG' : [9.1, 0.0235],
                    'GA' : [5.5, 0.0135],
                    'TC' : [8.6, 0.0229],
                    'CG' : [16.3, 0.0471],
                    'GC' : [8.0, 0.0171],
                    'GG' : [12.8, 0.0319],
                    'CC' : [9.3, 0.0232]
                        }
        for e in self.experiments:
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            fig_no = 0
            plot_no = 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                valleys = self.list_of_peaks[gene_name][e]['valleys']
                gene_end = self.genes[gene_name]['gene_length'] + self.five_prime_flank #as before slice_dataframe
                plot_no += 1
                fig.add_subplot(3, 2, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'])
                    if self.print_peaks == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['peaks_to_plot'], color='green')
                    if self.print_valleys == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['valleys_to_plot'], color='red')
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.4, color='orange')
                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e]['RT_to_plot'])
                    x_array = np.array(self.data[gene_name][e]['position'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                    if self.genes[gene_name]['RT'][e] is float():
                        plt.text(210-self.five_prime_flank,max(self.data[gene_name][e]['hits'])-150,'RT='+str(round(self.genes[gene_name]['RT'][e],3)))
                    else:
                        plt.text(210-self.five_prime_flank,max(self.data[gene_name][e]['hits'])-150,'RT=n/a')
                    plt.text(gene_end-200,max(self.data[gene_name][e]['hits'])-150,'delta G for: RNA/DNA_DNA/DNA')
                except KeyError:
                    plt.xlabel("NO READS")
                aaa = max(self.data[gene_name][e]['hits'])/10  #parameter to print energy in non-overlapping way
                inkr = 1
                print 'valleys: '+str(valleys)
                for v in valleys:
                    if v >= gene_end:
                        print 'for valley '+str(v)+' energy is:'
                        r_delta_H = 0
                        r_delta_S = 0
                        d_delta_H = 0
                        d_delta_S = 0
                        pair = str()
                        bbb = aaa * inkr
                        inkr += 1
                        for i in range(v-5,v+2):
                            print 'position: '+str(i)+' nucleotide: '+self.data[gene_name][e]['nucleotides'][i]
                            if len(pair) > 1:
                                pair = pair[1]
                            pair += self.data[gene_name][e]['nucleotides'][i]
                            if len(pair) == 2:
                                r_delta_H += rna_dna[pair][0]
                                r_delta_S += rna_dna[pair][1]
                                d_delta_H += dna_dna[pair][0]
                                d_delta_S += dna_dna[pair][1]
                        r_delta_G = r_delta_H * (1-(310/(r_delta_H/r_delta_S)))
                        d_delta_G = d_delta_H * (1-(310/(d_delta_H/d_delta_S)))
                        to_print = str(r_delta_G)+'_'+str(d_delta_G)
                        print 'r_delta_G: '+str(r_delta_G)
                        print 'd_delta_G: '+str(d_delta_G)
                        plt.text(v-self.five_prime_flank,bbb,to_print)
                if plot_no == 6:
                        fig_no += 1
                        plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png')
                        plt.clf()
                        plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    def termination_efficency(self):
        print '# Calculating energy for termination efficiency...'
        dna_dna = {	'AA' : [7.9, 0.0222],
                    'TT' : [7.9, 0.0222],
                    'AT' : [7.2, 0.0204],
                    'TA' : [7.2, 0.0213],
                    'CA' : [8.5, 0.0227],
                    'TG' : [8.5, 0.0227],
                    'GT' : [8.4, 0.0224],
                    'AC' : [8.4, 0.0224],
                    'CT' : [7.8, 0.021],
                    'AG' : [7.8, 0.021],
                    'GA' : [8.2, 0.0222],
                    'TC' : [8.2, 0.0222],
                    'CG' : [10.6, 0.0272],
                    'GC' : [9.8, 0.0244],
                    'GG' : [8.0, 0.0199],
                    'CC' : [8.0, 0.0199]
                        }
        rna_dna = { 'AA' : [7.8, 0.0219],
                    'TT' : [11.5, 0.0364],
                    'AT' : [8.3, 0.0239],
                    'TA' : [7.8, 0.0232],
                    'CA' : [9.0, 0.0261],
                    'TG' : [10.4, 0.0284],
                    'GT' : [7.8, 0.0216],
                    'AC' : [5.9, 0.0123],
                    'CT' : [7.0, 0.0197],
                    'AG' : [9.1, 0.0235],
                    'GA' : [5.5, 0.0135],
                    'TC' : [8.6, 0.0229],
                    'CG' : [16.3, 0.0471],
                    'GC' : [8.0, 0.0171],
                    'GG' : [12.8, 0.0319],
                    'CC' : [9.3, 0.0232]
                        }
        for e in self.experiments:
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            fig_no = 0
            plot_no = 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                gene_end = self.genes[gene_name]['gene_length'] + self.five_prime_flank #as before slice_dataframe
                plot_no += 1
                fig.add_subplot(3, 2, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'])
                    if self.print_peaks == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['peaks_to_plot'], color='green')
                    if self.print_valleys == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['valleys_to_plot'], color='red')
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.4, color='orange')
                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e]['RT_to_plot'])
                    x_array = np.array(self.data[gene_name][e]['position'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                    if self.genes[gene_name]['RT'][e] is float():
                        plt.text(210-self.five_prime_flank,max(self.data[gene_name][e]['hits'])-150,'RT='+str(round(self.genes[gene_name]['RT'][e],3)))
                    else:
                        plt.text(210-self.five_prime_flank,max(self.data[gene_name][e]['hits'])-150,'RT=n/a')
                    plt.text(gene_end-200,max(self.data[gene_name][e]['hits'])-150,'delta G for: RNA/DNA_DNA/DNA')
                except KeyError:
                    plt.xlabel("NO READS")
# calculating energy
                r_delta_H = 0
                r_delta_S = 0
                d_delta_H = 0
                d_delta_S = 0
                pair = str()
                print gene_name
                for i in range(gene_end +1,gene_end+21):
                    print 'position: '+str(i)+' nucleotide: '+self.data[gene_name][e]['nucleotides'][i]
                    if len(pair) > 1:
                        pair = pair[1]
                    pair += self.data[gene_name][e]['nucleotides'][i]
                    if len(pair) == 2:
                        r_delta_H += rna_dna[pair][0]
                        r_delta_S += rna_dna[pair][1]
                        d_delta_H += dna_dna[pair][0]
                        d_delta_S += dna_dna[pair][1]
                r_delta_G = r_delta_H * (1-(310/(r_delta_H/r_delta_S)))
                d_delta_G = d_delta_H * (1-(310/(d_delta_H/d_delta_S)))
                to_print = str(r_delta_G)+'_'+str(d_delta_G)
                print 'r_delta_G: '+str(r_delta_G)
                print 'd_delta_G: '+str(d_delta_G)
                plt.text(self.genes[gene_name]['gene_length'],(max(self.data[gene_name][e]['hits'])/2),to_print)
                if plot_no == 6:
                        fig_no += 1
                        plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png')
                        plt.clf()
                        plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    def slice_dataframe(self):
        five_prime_to_out = self.five_prime_flank - 50 # maybe changed - it simply gives only -50 flank independly what -r value was used for pyPileups
        for i in self.data:
            for e in self.data[i]:
                self.data[i][e] = self.data[i][e][five_prime_to_out::]
        # if self.list_of_peaks:
        #     for i in self.list_of_peaks:
        #         for e in self.list_of_peaks[i]:
        #             self.list_of_peaks[i][e]['peaks'] = [x - five_prime_to_out for x in self.list_of_peaks[i][e]['peaks']]
        #             self.list_of_peaks[i][e]['valleys'] = [x - five_prime_to_out for x in self.list_of_peaks[i][e]['valleys']]
        return True

### making output plots, one gene, different experiments per page
    def fig_gene_pp(self):
        print '# Plotting 1 gene per page (all experiments).'

        if len(self.experiments) <= 6:
            subplot_layout = [3, 2]
        elif len(self.experiments) > 6 and len(self.experiments) <= 9:
            subplot_layout = [3,3]
        elif len(self.experiments) > 9 and len(self.experiments) <= 12:
            subplot_layout = [4,3]
        elif len(self.experiments) > 12 and len(self.experiments) <= 16:
            subplot_layout = [4,4]
        else:
            print '# Unsupported layout - more than 16 different experiments'
            exit()

        for i_gene_id in self.genes_id_list:
            fig = plt.figure(figsize=(12, 9), dpi=200, facecolor='w', edgecolor='k')
            gene_name = self.id_to_names[i_gene_id]
            plot_no = 0
            for e in self.experiments:
                plot_no += 1
                fig.add_subplot(subplot_layout[0],subplot_layout[1],plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'])
                    if self.print_peaks == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['peaks_to_plot'])
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.4, color='orange')
                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e]['RT_to_plot'])
                    x_array = np.array(self.data[gene_name][e]['position'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                except KeyError:
                    plt.xlabel("NO READS")

                if plot_no == len(self.experiments):
                    plot_no = 0
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'.png', dpi=200)
                    plt.clf()
        return True
###function to plot ratio between sample with -a parameter and -b parameter
    def fig_ratio(self, a, b):
        print '# Plotting 1 gene per page. Log2 ratio for '+a+' divided by '+b+' (all experiments).'
        new_exp_list = self.group_experiments(a,b)
        if len(new_exp_list) <= 6:
            subplot_layout = [3, 2]
        elif len(new_exp_list) > 6 and len(new_exp_list) <= 9:
            subplot_layout = [3,3]
        elif len(new_exp_list) > 9 and len(new_exp_list) <= 12:
            subplot_layout = [4,3]
        elif len(new_exp_list) > 12 and len(new_exp_list) <= 16:
            subplot_layout = [4,4]
        else:
            print '# fig_ratio: Unsupported layout - more than 16 different experiments'
            exit()

        for i_gene_id in self.genes_id_list:
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            gene_name = self.id_to_names[i_gene_id]
            plot_no = 0
            for e in new_exp_list:
                plot_no += 1
                fig.add_subplot(subplot_layout[0],subplot_layout[1],plot_no)
                plt.tight_layout()
                plt.title('log2 '+e[0]+'/'+e[1]+' ratio')
                plt.ylabel("no. of reads")
                try:
                    self.data[gene_name][e[0]]['hits'] = np.log2(self.data[gene_name][e[0]]['hits']/self.data[gene_name][e[1]]['hits']) ###getting ratio
                    plt.plot(self.data[gene_name][e[0]]['position'], self.data[gene_name][e[0]]['hits'])
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.4, color='orange')
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                except KeyError:
                    plt.xlabel("NO READS")
                if plot_no == len(new_exp_list):
                    plot_no = 0
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'.png')
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
            gene = re.search(r'(.*)_'+a, e).group(1)
            second = [d for d in b_experiments if re.search(gene,d)]
            if len(second) > 1:
                print 'To many elements to divide for '+gene
                exit()
            couple = [e, second[0]]
            ratio_exp_list.append(couple)
        print 'Found pairs:'
        print ratio_exp_list
        return ratio_exp_list

# print plots, gene after gene,
    def fig_gene_after_gene(self):
        print '# Plotting gene after gene.'
        for e in self.experiments:
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            fig_no = 0
            plot_no = 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                plot_no += 1
                fig.add_subplot(3, 2, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'])
                    if self.print_peaks == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['peaks_to_plot'], color='green')
                    if self.print_valleys == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['valleys_to_plot'], color='red')
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.4, color='orange')
                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e]['RT_to_plot'])
                    x_array = np.array(self.data[gene_name][e]['position'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                except KeyError:
                    plt.xlabel("NO READS")
                if plot_no == 6:
                    fig_no += 1
                    plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png', dpi=200)
                    plt.clf()
                    plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no+1)+'.png', dpi=200)
                plt.clf()
        return True

# print plots, gene after gene,
    def fig_boxes(self,abox,bbox):
        print '# Reading A and A box files'
        abox_dict = dict()
        abox_len = 16
        for line in abox:
            try:
                if line[0] == "#": continue
                line_elements = line.strip().split('\t')
                gene_name, start = line_elements[0],int(line_elements[1])
                if gene_name not in abox_dict:
                    abox_dict[gene_name] = start
            except IndexError:
                sys.stderr.write("\nIndexError at line:\n%s\n" % line)
                pass

        bbox_dict = dict()
        bbox_len = 9
        for line in bbox:
            try:
                if line[0] == "#": continue
                line_elements = line.strip().split('\t')
                gene_name, start = line_elements[0],int(line_elements[1])
                if gene_name not in bbox_dict:
                    bbox_dict[gene_name] = start
            except IndexError:
                sys.stderr.write("\nIndexError at line:\n%s\n" % line)
                pass

        print '# Plotting gene after gene, marking boxes.'
        for e in self.experiments:
            fig = plt.figure(figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
            fig_no = 0
            plot_no = 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                plot_no += 1
                fig.add_subplot(3, 2, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'])
                    if self.print_peaks == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['peaks_to_plot'], color='green')
                    if self.print_valleys == True:
                        plt.plot(self.data[gene_name][e]['position'], self.data[gene_name][e]['valleys_to_plot'], color='red')
            #plot exons background
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.2, color='orange')
            #plot A and B box
                    try:
                        plt.axvspan(abox_dict[gene_name], abox_dict[gene_name]+abox_len, alpha=0.6, color='red')
                    except:
                        pass
                    try:
                        plt.axvspan(bbox_dict[gene_name], bbox_dict[gene_name]+bbox_len, alpha=0.6, color='green')
                    except:
                        pass
                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e]['RT_to_plot'])
                    x_array = np.array(self.data[gene_name][e]['position'])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                except KeyError:
                    plt.xlabel("NO READS")
                if plot_no == 6:
                    fig_no += 1
                    plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png')
                    plt.clf()
                    plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

# print histogram for 3` end including nucleotides
    def fig_nucleotide_resolution(self):
        print '# Plotting 3` end.'
        fig = plt.figure(figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
        if len(self.experiments) > 1:
            print 'More than one experiment. Preprocess concat file to get only one experiment:' \
                  "cat your.concat | awk '$7 == \"experiment\" {print $0}' > experiment.concat"
            exit()
        else:
            e = self.experiments[0]
        for tRNA_group in self.rt:
            fig_no = 0
            plot_no = 0
            rt_sorted = collections.OrderedDict(sorted(self.rt[tRNA_group].items(), reverse=True, key=lambda t: t[1]))
            for i_gene_id in rt_sorted:
                gene_name = self.id_to_names[i_gene_id]
                plot_no += 1
                ax = fig.add_subplot(5, 1, plot_no)
                fig.tight_layout()
                plt.title(tRNA_group)
                try:
                    # self.data[gene_name][e] = self.data[gene_name][e][(-10-self.three_prime_flank):(120-self.three_prime_flank):] # plot from last 130 nt
                    # print self.data[gene_name][e]
                    # bar = ax.bar(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'], width=0.5)
                    # for i in range(0,10):
                    #     bar[i].set_color('green')
                    # for i in range(10,len(bar)):
                    #     bar[i].set_color('grey')
                    # plt.xticks(list(self.data[gene_name][e]['position']), list(self.data[gene_name][e]['nucleotides']), fontsize=8)
                    uno_to_plot = self.data[gene_name][e][(-10-self.three_prime_flank):(120-self.three_prime_flank):] # plot from last 130 nt
                    ax.set_ylabel("no. of reads")
                    ax.set_xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                    bar = ax.bar(uno_to_plot['position'], uno_to_plot['hits'], width=0.5)
                    axes = plt.gca()
                    axes.set_ylim([0,2000])
                    for i in range(0,10):
                        bar[i].set_color('green')
                    for i in range(10,len(bar)):
                        bar[i].set_color('grey')
                    plt.xticks(list(uno_to_plot['position']), list(uno_to_plot['nucleotides']), fontsize=8)

                    plot_no += 1
                    ax = fig.add_subplot(5, 1, plot_no)
                    ax.set_ylabel("no. of reads")
                    ax.set_xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                    duo_to_plot = self.data[gene_name][e][(120-self.three_prime_flank)::] # plot from last 130 nt
                    bar = ax.bar(duo_to_plot['position'], duo_to_plot['hits'], width=0.5, color='grey')
                    axes = plt.gca()
                    axes.set_ylim([0,2000])
                    # axes.set_ylim([0,max(self.data[gene_name][e]['hits'])])
                    plt.xticks(list(duo_to_plot['position']), list(duo_to_plot['nucleotides']), fontsize=8)

                except KeyError:
                    plt.text(0.5,0.5,"NO READS")
                if plot_no == 5:
                    fig_no += 1
                    plt.savefig(self.prefix+'nuc'+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_'+tRNA_group+'_fig_'+str(fig_no)+'.png')
                    plt.clf()
                    plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+'nuc'+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_'+tRNA_group+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    def fig_nucleotide_gene(self):
        print '# Plotting gene in nucleotide resolution.'
        fig = plt.figure(figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
        if len(self.experiments) > 1:
            print 'More than one experiment. Preprocess concat file to get only one experiment:' \
                  "cat your.concat | awk '$7 == \"experiment\" {print $0}' > experiment.concat"
            exit()
        else:
            e = self.experiments[0]
            fig_no, plot_no = 0, 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                plot_no += 1
                ax = fig.add_subplot(5, 1, plot_no)
                fig.tight_layout()
                plt.title(e)
                ax.set_ylabel("no. of reads")
                ax.set_xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                try:
                    self.data[gene_name][e] = self.data[gene_name][e][self.three_prime_flank:-self.three_prime_flank:] # plot only gene
                    bar = ax.bar(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'], width=0.5)
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.2, color='orange')
                    plt.xticks(list(self.data[gene_name][e]['position']), list(self.data[gene_name][e]['nucleotides']), fontsize=8)
                except KeyError:
                    plt.text(0.5,0.5,"NO READS")
                if plot_no == 5:
                    fig_no += 1
                    plt.savefig(self.prefix+'nuc_gene'+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_'+'_fig_'+str(fig_no)+'.png')
                    plt.clf()
                    plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+'nuc_gene'+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_'+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    def fig_energy(self, step):
        print 'Plotting 3` end and calculating energy for plots...'
        dna_dna = {	'AA' : [7.9, 0.0222],
                    'TT' : [7.9, 0.0222],
                    'AT' : [7.2, 0.0204],
                    'TA' : [7.2, 0.0213],
                    'CA' : [8.5, 0.0227],
                    'TG' : [8.5, 0.0227],
                    'GT' : [8.4, 0.0224],
                    'AC' : [8.4, 0.0224],
                    'CT' : [7.8, 0.021],
                    'AG' : [7.8, 0.021],
                    'GA' : [8.2, 0.0222],
                    'TC' : [8.2, 0.0222],
                    'CG' : [10.6, 0.0272],
                    'GC' : [9.8, 0.0244],
                    'GG' : [8.0, 0.0199],
                    'CC' : [8.0, 0.0199]
                        }
        rna_dna = { 'AA' : [7.8, 0.0219],
                    'TT' : [11.5, 0.0364],
                    'AT' : [8.3, 0.0239],
                    'TA' : [7.8, 0.0232],
                    'CA' : [9.0, 0.0261],
                    'TG' : [10.4, 0.0284],
                    'GT' : [7.8, 0.0216],
                    'AC' : [5.9, 0.0123],
                    'CT' : [7.0, 0.0197],
                    'AG' : [9.1, 0.0235],
                    'GA' : [5.5, 0.0135],
                    'TC' : [8.6, 0.0229],
                    'CG' : [16.3, 0.0471],
                    'GC' : [8.0, 0.0171],
                    'GG' : [12.8, 0.0319],
                    'CC' : [9.3, 0.0232]
                        }
        fig = plt.figure(figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
        if len(self.experiments) > 1:
            print 'More than one experiment. Preprocess concat file to get only one experiment:' \
                  "cat your.concat | awk '$7 == \"experiment\" {print $0}' > experiment.concat"
            exit()
        else:
            e = self.experiments[0]
        for tRNA_group in self.rt:
            fig_no = 0
            plot_no = 0
            rt_sorted = collections.OrderedDict(sorted(self.rt[tRNA_group].items(), reverse=True, key=lambda t: t[1]))
            for i_gene_id in rt_sorted:
                gene_name = self.id_to_names[i_gene_id]
                plot_no += 1
                fig.add_subplot(5, 1, plot_no)
                fig.tight_layout()
                plt.title(tRNA_group)
                plt.ylabel("no. of reads")
                plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                plt.text(str(self.genes[gene_name]['RT'][e]), horizontalalignment='right', verticalalignment='top')
                try:
                    self.data[gene_name][e] = self.data[gene_name][e][(-10-self.three_prime_flank):(140-self.three_prime_flank):] # plot from last 150 nt
                    bar = plt.bar(self.data[gene_name][e]['position'], self.data[gene_name][e]['hits'], width=0.6)
                    for i in range(0,10):
                        bar[i].set_color('green')
                    for i in range(10,len(bar)):
                        bar[i].set_color('grey')
                    plt.xticks(list(self.data[gene_name][e]['position']), list(self.data[gene_name][e]['nucleotides']), fontsize=8)
                    ax2 = plt.twinx()
                    ax2.set_ylabel("energy")

                ## plot energy every five nucleotides
                    r_delta_H = 0
                    r_delta_S = 0
                    d_delta_H = 0
                    d_delta_S = 0
                    window = step
                    pair = str()
                    list_for_G = list()
                    for i in range(min(self.data[gene_name][e].index),max(self.data[gene_name][e].index)+1):
                        if len(pair) > 1:
                            pair = pair[1]
                        pair += self.data[gene_name][e]['nucleotides'][i]
                        if len(pair) == 2:
                            r_delta_H += rna_dna[pair][0]
                            r_delta_S += rna_dna[pair][1]
                            d_delta_H += dna_dna[pair][0]
                            d_delta_S += dna_dna[pair][1]
                            list_for_G.append(i)
                        window -= 1
                        if window == 0:
                            r_delta_G = r_delta_H * (1-(310/(r_delta_H/r_delta_S)))
                            d_delta_G = d_delta_H * (1-(310/(d_delta_H/d_delta_S)))
                            for a in list_for_G:
                                self.data[gene_name][e].loc[a,'r_delta_G'] = r_delta_G
                                self.data[gene_name][e].loc[a,'d_delta_G'] = d_delta_G

                            list_for_G = []
                            r_delta_G = 0
                            d_delta_G = 0
                            r_delta_H = 0
                            r_delta_S = 0
                            d_delta_H = 0
                            d_delta_S = 0
                            window = step
                    plt.plot(list(self.data[gene_name][e]['position']), list(self.data[gene_name][e]['r_delta_G']))
                    plt.plot(list(self.data[gene_name][e]['position']), list(self.data[gene_name][e]['d_delta_G']))

                except KeyError:
                    plt.text(10, 10, "NO READS")
                if plot_no == 5:
                    fig_no += 1
                    plt.savefig(self.prefix+'nuc_energy'+'_w'+str(step)+'_'+tRNA_group+'_fig_'+str(fig_no)+'.png')
                    plt.clf()
                    plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+'nuc_energy'+'_w'+str(step)+'_'+tRNA_group+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
        return True

    def test_print(self, what):
        pd.set_option('display.max_rows', 5000)
        print what
        pd.reset_option('display.max_rows')
        return True