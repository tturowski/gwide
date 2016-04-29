#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "2.0"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

import numpy as np
import sys, collections, re
from pypeaks import Data
from pyCRAC.Parsers import GTF2
import matplotlib.pyplot as plt
import pandas as pd

class mRNAFromConcat():
    def __init__(self, gtf_file, five_prime_flank, three_prime_flank, hits_threshold, lookahead, prefix, npM):
        self.gtf_file           =   str(gtf_file)
        self.gtf                =   GTF2.Parse_GTF()
        self.gtf.read_GTF(self.gtf_file)
        self.genes              =   dict()
        self.data               =   dict()
        self.id_to_names        =   dict()
        # self.rt                 =   dict()                  # designated to work with one experiment only
        self.genes_name_list    =   list()
        self.genes_id_list      =   list()
        self.five_prime_flank   =   five_prime_flank
        self.three_prime_flank  =   three_prime_flank
        # self.five_prime_to_print=   50
        self.longest_gene       =   0
        self.hits_threshold     =   hits_threshold
        self.lookahead          =   lookahead
        self.prefix             =   str(prefix)
        # self.print_peaks        =   print_peaks
        # self.print_valleys      =   print_valleys
        # self.readthrough_start  =   readthrough_start
        self.experiments        =   list()
        self.npM         =   npM  # -n option allows for using normalized dataset (reads p[er Million)
        self.list_of_peaks = dict()
        self.dna_dna = {	'AA' : [7.9, 0.0222],
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
        self.rna_dna = {    'AA' : [7.8, 0.0219],
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

    def read_csv(self, concat_file,use='hits'):
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
                # 'introns'       :   self.get_introns(gene_name),   # [[len,len...][(start,stop),etc]] <-[[list_with_lengths],[list_with_start and stop]]
                # 'RT'            :   dict(),
                                    }
            # self.genes[gene_name]['exons'] = self.get_exons(gene_name)
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
                # print gene_csv[gene_csv.exp_name == e] # to check
                self.data[gene_name][e] = gene_csv[gene_csv.exp_name == e][use]
                if self.npM == True:
                    self.data[gene_name][(e+"_nrpm")] = gene_csv[gene_csv.exp_name == e]['n_'+use]    #to work with data normalized reads per M
            self.data[gene_name] = self.data[gene_name].fillna(0)
            # print self.data[gene_name] #line to check
        self.genes_id_list.sort()
        return True

    def bind(self, exp_to_use,window):
        """find binding fragment using deletions from CRAC data"""
        for gene_name in self.genes_name_list:
            gene_data = self.data[gene_name]
            max_position = gene_data[exp_to_use].idxmax()
            motif = ''.join(list(gene_data[max_position-window:max_position+window]['nucleotide']))
            if len(motif) >= 8:
                print '>'+gene_name
                print motif + '\n'
                # print gene_data[max_position-window:max_position+window]

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
                print i
                print e
                print list(self.data[i][e])
                hist = Data(list(self.data[i].index), list(self.data[i][e]), smoothness=1, default_smooth=False)
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
        print self.list_of_peaks
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
                    self.genes[gene_name]['total']      =   total
                    self.genes[gene_name]['total_av']   =   total_av
                    self.genes[gene_name]['total_med']  =   total_med
                    self.genes[gene_name]['total_std']  =   total_SD
                    self.genes[gene_name]['a']          =   a
                    self.genes[gene_name]['b']          =   b
                    self.genes[gene_name]['i']          =   i
                    self.genes[gene_name]['e']          =   e
                    self.genes[gene_name]['f']          =   f
                    self.genes[gene_name]['d']      =   d
                    self.genes[gene_name]['d_av']   =   d_av
                    self.genes[gene_name]['d_med']  =   d_med
                    self.genes[gene_name]['d_std']  =   d_SD
                    self.genes[gene_name]['g']      =   g
                    self.genes[gene_name]['g_av']   =   g_av
                    self.genes[gene_name]['g_med']  =   g_med
                    self.genes[gene_name]['g_std']  =   g_SD

        ### below part of function is needed to sort tRNA according to readthrough
                tRNA_group = self.genes[gene_name]['gene_id'][0:2]
                gene_id = self.genes[gene_name]['gene_id']
                if tRNA_group not in self.rt:
                    self.rt[tRNA_group] = dict()
                self.rt[tRNA_group][gene_id] = RT
        return True

# making output; if -p option then text file, else to standard output
    def make_text_file(self, filename, details=False, print_dG=False, ntotal=False, nmax=False):
        print '# Making text file...'
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
            if print_dG == True:
                output_dict[i].append(str(self.genes[i]['r_delta_G']))
                output_dict[i].append(str(self.genes[i]['d_delta_G']))
            for e in self.experiments:
                output_dict[i].append(str(self.genes[i]['RT'][e]))
            if details==True:
                output_dict[i].append(str(self.genes[i]['total']))
                output_dict[i].append(str(self.genes[i]['total_av']))
                output_dict[i].append(str(self.genes[i]['total_med']))
                output_dict[i].append(str(self.genes[i]['total_std']))
                output_dict[i].append(str(self.genes[i]['a']))
                output_dict[i].append(str(self.genes[i]['b']))
                output_dict[i].append(str(self.genes[i]['i']))
                output_dict[i].append(str(self.genes[i]['e']))
                output_dict[i].append(str(self.genes[i]['f']))
                output_dict[i].append(str(self.genes[i]['d']))
                output_dict[i].append(str(self.genes[i]['d_av']))
                output_dict[i].append(str(self.genes[i]['d_med']))
                output_dict[i].append(str(self.genes[i]['d_std']))
                output_dict[i].append(str(self.genes[i]['g']))
                output_dict[i].append(str(self.genes[i]['g_av']))
                output_dict[i].append(str(self.genes[i]['g_med']))
                output_dict[i].append(str(self.genes[i]['g_std']))
        # experiments = '\t'.join(self.experiments)
        filename.write("# analyse_tRNA_pileups output file using gtf file:"+"\n")
        filename.write("# "+self.gtf_file+"\n")
        filename.write("# 5` flank: "+str(self.five_prime_flank)+"\n")
        filename.write("# 3` flank: "+str(self.three_prime_flank)+"\n")
        # if self.normalized == True:
        #     filename.write("# calculations performed on data normalized to reads per Million \n")
        # else:
        #     filename.write("# calculations performed on non-normalized data \n")
        filename.write("# readthrough calculated: last nucleotide of gene + "+str(self.readthrough_start)+" nt\n")
        filename.write("# threshold of hits: "+str(self.hits_threshold)+"\n")
        filename.write("# lookahead option for pypeaks: "+str(self.lookahead)+"\n")
        if nmax == True and ntotal == False:
            filename.write("# Calculations performed on normalized data (max = 1)"+"\n")
        elif ntotal == True:
            filename.write("# Calculations performed on normalized data (sum = 1)"+"\n")
        else:
            filename.write("# Calculations performed on non-normalized data "+"\n")
        header_list = ["# name","ID","length","int_len"]
        if print_dG == True:
            header_list += ["dG_RNA/DNA","dG_DNA/DNA"]
        header_list += self.experiments
        if details == True:
            header_list += ['total','total_av','total_med','total_std','a','b','i','e','f','d','d_av','d_med','d_std','g','g_av','g_med','g_std']
        filename.write('\t'.join(header_list)+'\n')
        for g in self.genes_id_list:
            line = self.id_to_names[g]
            filename.write('\t'.join(output_dict[line]) + '\n')
        if filename != sys.stdout:
            filename.close()
        return True

    def slice_dataframe(self):
        five_prime_to_out = self.five_prime_flank - self.five_prime_to_print # maybe changed - it simply gives only -50 flank independly what -r value was used for pyPileups
        for i in self.data:
            for e in self.data[i]:
                self.data[i][e] = self.data[i][e][five_prime_to_out::]
        # if self.list_of_peaks:
        #     for i in self.list_of_peaks:
        #         for e in self.list_of_peaks[i]:
        #             self.list_of_peaks[i][e]['peaks'] = [x - five_prime_to_out for x in self.list_of_peaks[i][e]['peaks']]
        #             self.list_of_peaks[i][e]['valleys'] = [x - five_prime_to_out for x in self.list_of_peaks[i][e]['valleys']]
        return True