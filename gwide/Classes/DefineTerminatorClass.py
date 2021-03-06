#!/usr/bin/env python
import numpy as np
import sys, collections, re
from pypeaks import Data
from pyCRAC.Parsers import GTF2
import matplotlib.pyplot as plt
import pandas as pd

class DefineTerminatorClass():
    def __init__(self, gtf_file, five_prime_flank, three_prime_flank, hits_threshold, lookahead, prefix, print_valleys, print_peaks, readthrough_start, normalized):
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
        self.five_prime_to_print=   50
        self.longest_gene       =   0
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

    def read_csv(self, concat_file, null_substitution=False):
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
            self.data[gene_name]['nucleotide'] = gene_csv['nucleotide'][:(self.five_prime_flank + gene_length + self.three_prime_flank):]
            for e in self.experiments:
                self.data[gene_name][e] = gene_csv[gene_csv.exp_name == e]['hits']
                if self.normalized == True:
                    self.data[gene_name][(e+"_nrpm")] = gene_csv[gene_csv.exp_name == e]['n_hits']    #to work with data normalized reads per M
            self.data[gene_name] = self.data[gene_name].fillna(0)
            if null_substitution == True:
                self.data[gene_name] = self.data[gene_name].replace(0,1)
            # print self.data[gene_name][230:330]
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
                    hist.get_peaks(method="slope", peak_amp_thresh=0.00005, valley_thresh=0.00003, intervals=None,
                                   lookahead=self.lookahead, avg_interval=100)
                    self.list_of_peaks[i][e]['peaks'] = sorted(np.array(hist.peaks['peaks'][0]).tolist())
                    self.list_of_peaks[i][e]['valleys'] = sorted(np.array(hist.peaks['valleys'][0]).tolist())
                except ValueError:
                    pass
                    self.list_of_peaks[i][e]['peaks'] = []
                    self.list_of_peaks[i][e]['valleys'] = []
        # print self.list_of_peaks
        return True

    def calculate(self, details=False, ntotal=False, nmax=False):
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
                    histogram   =   list(self.data[gene_name][exp])
                    if max(histogram) >= self.hits_threshold:
                #normalization options
                        if nmax == True:
                            exp = exp_old+'_nmax'
                            self.data[gene_name][exp] = self.data[gene_name][exp_old]/self.data[gene_name][exp_old].max()
                        if ntotal == True:
                            exp = exp_old+'_ntotal'
                            self.data[gene_name][exp] = self.data[gene_name][exp_old]/self.data[gene_name][exp_old].sum()
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

    def calculate_dG(self):
        print '# Calculating delta G energy for RNA/DNA and DNA/DNA for 20 nt after 3` end'
        for gene_name in self.genes_name_list:
            gene_end = self.genes[gene_name]['gene_length'] + self.five_prime_flank #as before slice_dataframe
        # calculating energy
            r_delta_H = 0
            r_delta_S = 0
            d_delta_H = 0
            d_delta_S = 0
            pair = str()
            for i in range(gene_end+1,gene_end+21):
                if len(pair) > 1:
                    pair = pair[1]
                pair += self.data[gene_name]['nucleotide'][i]
                if len(pair) == 2:
                    r_delta_H += self.rna_dna[pair][0]
                    r_delta_S += self.rna_dna[pair][1]
                    d_delta_H += self.dna_dna[pair][0]
                    d_delta_S += self.dna_dna[pair][1]
            r_delta_G = r_delta_H * (1-(310/(r_delta_H/r_delta_S)))
            d_delta_G = d_delta_H * (1-(310/(d_delta_H/d_delta_S)))
            # to_print = str(r_delta_G)+'_'+str(d_delta_G)
            self.genes[gene_name]['r_delta_G'] = r_delta_G
            self.genes[gene_name]['d_delta_G'] = d_delta_G
        return True

    #calculating T richness for each nucleotide together with 8 nucleotides before
    def calculate_T_richness(self):
        print "# Calculating T richness..."
        #to avoid calculations for fist 200 nt
        five_prime_to_out = self.five_prime_flank - self.five_prime_to_print
        for gene_name in self.genes_name_list:
            for index, nucleotide in self.data[gene_name]['nucleotide'][five_prime_to_out::].iteritems():
                #getting 9 nt including present and picking up only Ts
                temp_data = self.data[gene_name][index-9:index:]
                T_richness = len(temp_data[temp_data.nucleotide == "T"])
                #saving no of T per 9 nt as T richness
                self.data[gene_name].loc[index, "Trichness"] = T_richness

    def find_terminators(self):
        def flat_list(list_to_flat):
            new_list = [item for sublist in list_to_flat for item in sublist]
            return new_list

        #extending found terminators according to their length
        def extend(position, no_to_extend):
            extended = [position + no for no in range(0,no_to_extend)]
            return extended

        #finding all T streches, saves start, stop and previous and next nucleotides in dataframe
        def find_Tstretch(sequence, gene_name):
            #getting positions and lengths of terminators
            terminators_dict = dict()
            busy_list = list()
            data = pd.DataFrame(columns=["type", "type_numeric", "start", "stop", "bp_previous", "bp_next", "type_previous", "type_previous_distance", "type_next", "type_next_distance"])
            for no in range(20,1,-1): #not counting T1
                pattern = no * "T"
                type = "T"+str(no)
                terminators_list = [term.start() for term in re.finditer(pattern, sequence)]
                #remove identified previously
                terminators_list = [position for position in terminators_list if position not in busy_list]
                #extend and add to busy list
                terminators_extended = [extend(position=position, no_to_extend=no) for position in terminators_list]
                for position, start in enumerate(terminators_list):
                    stop = start+no-1 # -1 to make it human readable, but positions starts from 0
                    if stop+2 > len(sequence):
                        row = pd.DataFrame([[gene_name, type, no, start, stop, sequence[start-1], "out_of_index"]], columns=["gene_name", "type", "type_numeric", "start", "stop", "bp_previous", "bp_next"])
                    elif start == 0:
                        row = pd.DataFrame([[gene_name, type, no, start, stop, "out_of_index", sequence[stop+1]]], columns=["gene_name", "type", "type_numeric", "start", "stop", "bp_previous", "bp_next"])
                    else:
                        row = pd.DataFrame([[gene_name, type, no, start, stop, sequence[start-1], sequence[stop+1]]], columns=["gene_name", "type", "type_numeric", "start", "stop", "bp_previous", "bp_next"])
                    data = data.append(row, ignore_index=True)
                busy_list += flat_list(terminators_extended)
                terminators_dict["T"+str(no)] = terminators_list
            return data

        #analyzing previous dataframe to find neighbours of terminators
        def define_neighborhood(terminators_df):
            for index, row in terminators_df.iterrows():
                #find previous terminator
                previous_Tstretch = terminators_df[terminators_df.stop < row.start]['stop']
                if len(previous_Tstretch) > 0:
                    previous_Tstretch_index = previous_Tstretch.argmax()
                    previous_Tstretch_row = terminators_df.loc[previous_Tstretch_index]
                    terminators_df.loc[index, "type_previous"] = previous_Tstretch_row["type_numeric"]
                    terminators_df.loc[index, "type_previous_distance"] = row.start - previous_Tstretch_row["stop"]

                #find next terminator
                next_Tstretch = terminators_df[terminators_df.start > row.stop]['start']
                if len(next_Tstretch) > 0:
                    next_Tstretch_index = next_Tstretch.argmin()
                    next_Tstretch_row = terminators_df.loc[next_Tstretch_index]
                    terminators_df.loc[index, "type_next"] = next_Tstretch_row["type_numeric"]
                    terminators_df.loc[index, "type_next_distance"] = next_Tstretch_row["start"] - row.stop
            return terminators_df

        #calculate efficiency for each terminator
        def termination_efficiency(data, exp, terminators_df, length_to_add):
            terminators_df = terminators_df[terminators_df.type_numeric > 2] #calculate only for T3 and longer
            for index, row in terminators_df.iterrows():
                terminator_position = int(length_to_add + row.start + round(row.type_numeric/2, 0))
                for i in [30,40,50,60]:
                    previous_reads = data[terminator_position-i:terminator_position].mean()[exp]
                    next_reads = data[terminator_position:terminator_position+i].mean()[exp]
                    terminators_df.loc[index, "range_"+str(i)] = np.float64(next_reads)/previous_reads # np.float64 to avoid problems with 0
            return terminators_df

        print "# Finding terminators..."
        for e in self.experiments:
            terminators_all = pd.DataFrame()
            a = 0
            for gene_name in self.genes_name_list:
                a += 1
                print "# Calculating "+gene_name+" ("+str(a)+" of "+str(len(self.genes_name_list))+")"
                #making 3` end sequence
                sequence = "".join(self.data[gene_name]['nucleotide'][-260::].tolist()) #last 10 nt of gene and 3` end
                # print sequence
                terminators_df = find_Tstretch(sequence=sequence, gene_name=gene_name)
                terminators_df = define_neighborhood(terminators_df)
                length_to_add = self.five_prime_flank + self.genes[gene_name]['gene_length'] - 10
                # length_to_add = self.data[gene_name]['position'].max()-10
                terminators_df = termination_efficiency(data=self.data[gene_name], exp=e, terminators_df=terminators_df, length_to_add=length_to_add)
                terminators_all = terminators_all.append(terminators_df)
            terminators_all.to_csv("terminators_"+e+".table", sep='\t', index=False)

# making output; if -p option then text file, else to standard output

    def first_terminator(self):
        def flat_list(list_to_flat):
            new_list = [item for sublist in list_to_flat for item in sublist]
            return new_list

        #extending found terminators according to their length
        def extend(position, no_to_extend):
            extended = [position + no for no in range(0,no_to_extend)]
            return extended

        #finding all T streches, saves start, stop and previous and next nucleotides in dataframe
        def find_Tstretch(sequence, gene_name):
            #getting positions and lengths of terminators
            terminators_dict = dict()
            busy_list = list()
            data = pd.DataFrame(columns=["type", "type_numeric", "start", "stop", "bp_previous", "bp_next", "type_previous", "type_previous_distance", "type_next", "type_next_distance"])
            for no in range(20,1,-1): #not counting T1
                pattern = no * "T"
                type = "T"+str(no)
                terminators_list = [term.start() for term in re.finditer(pattern, sequence)]
                #remove identified previously
                terminators_list = [position for position in terminators_list if position not in busy_list]
                #extend and add to busy list
                terminators_extended = [extend(position=position, no_to_extend=no) for position in terminators_list]
                for position, start in enumerate(terminators_list):
                    stop = start+no-1 # -1 to make it human readable, but positions starts from 0
                    if stop+2 > len(sequence):
                        row = pd.DataFrame([[gene_name, type, no, start, stop, sequence[start-1], "out_of_index"]], columns=["gene_name", "type", "type_numeric", "start", "stop", "bp_previous", "bp_next"])
                    elif start == 0:
                        row = pd.DataFrame([[gene_name, type, no, start, stop, "out_of_index", sequence[stop+1]]], columns=["gene_name", "type", "type_numeric", "start", "stop", "bp_previous", "bp_next"])
                    else:
                        row = pd.DataFrame([[gene_name, type, no, start, stop, sequence[start-1], sequence[stop+1]]], columns=["gene_name", "type", "type_numeric", "start", "stop", "bp_previous", "bp_next"])
                    data = data.append(row, ignore_index=True)
                busy_list += flat_list(terminators_extended)
                terminators_dict["T"+str(no)] = terminators_list
            return data

        #analyzing previous dataframe to find neighbours of terminators
        def define_neighborhood(terminators_df):
            for index, row in terminators_df.iterrows():
                #find previous terminator
                previous_Tstretch = terminators_df[terminators_df.stop < row.start]['stop']
                if len(previous_Tstretch) > 0:
                    previous_Tstretch_index = previous_Tstretch.argmax()
                    previous_Tstretch_row = terminators_df.loc[previous_Tstretch_index]
                    terminators_df.loc[index, "type_previous"] = previous_Tstretch_row["type_numeric"]
                    terminators_df.loc[index, "type_previous_distance"] = row.start - previous_Tstretch_row["stop"]

                #find next terminator
                next_Tstretch = terminators_df[terminators_df.start > row.stop]['start']
                if len(next_Tstretch) > 0:
                    next_Tstretch_index = next_Tstretch.argmin()
                    next_Tstretch_row = terminators_df.loc[next_Tstretch_index]
                    terminators_df.loc[index, "type_next"] = next_Tstretch_row["type_numeric"]
                    terminators_df.loc[index, "type_next_distance"] = next_Tstretch_row["start"] - row.stop
            return terminators_df

        #calculate efficiency for each terminator
        def termination_efficiency(data, exp, gene_name, terminators_df, length_to_add):
            terminators_df = terminators_df[terminators_df.type_numeric > 2] #calculate only for T3 and longer
            terminators_df = terminators_df[terminators_df.start == terminators_df['start'].min()] #calculate only for first terminator
            for index, row in terminators_df.iterrows():
                terminator_position = int(length_to_add + row.start + round(row.type_numeric/2, 0))
                for i in [30,40,50,60]:
                    previous_reads = data[gene_name][terminator_position-i:terminator_position].mean()[exp]
                    next_reads = data[gene_name][terminator_position:terminator_position+i].mean()[exp]
                    terminators_df.loc[index, "range_"+str(i)] = np.float64(next_reads)/previous_reads # np.float64 to avoid problems with 0
            return terminators_df

        print "# Finding terminators..."
        for e in self.experiments:
            terminators_all = pd.DataFrame()
            a = 0
            for gene_name in self.genes_name_list:
                a += 1
                print "# Calculating "+gene_name+" ("+str(a)+" of "+str(len(self.genes_name_list))+")"
                #making 3` end sequence
                sequence = "".join(self.data[gene_name]['nucleotide'][-260::].tolist()) #last 10 nt of gene and 3` end
                # print sequence
                terminators_df = find_Tstretch(sequence=sequence, gene_name=gene_name)
                terminators_df = define_neighborhood(terminators_df)
                length_to_add = self.five_prime_flank + self.genes[gene_name]['gene_length'] - 10
                # length_to_add = self.data[gene_name]['position'].max()-10
                terminators_df = termination_efficiency(data=self.data, exp=e, gene_name=gene_name, terminators_df=terminators_df, length_to_add=length_to_add)
                terminators_all = terminators_all.append(terminators_df)
            terminators_all.to_csv("terminators_"+e+".table", sep='\t', index=False)

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

    def fig_Trichness(self):
        print '# Plotting and T richness for 9 nucleotides...'
        for e in self.experiments:
            fig = plt.figure(figsize=(12, 9), dpi=200, facecolor='w', edgecolor='k')
            fig_no = 0
            plot_no = 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                gene_length =   self.genes[gene_name]['gene_length']
                RT_begin    =   self.five_prime_flank + gene_length + self.readthrough_start
                plot_no += 1
                fig.add_subplot(3, 2, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                try:
                    plt.plot(self.data[gene_name]['position'], self.data[gene_name][e], color='black')
                    if self.print_peaks == True:
                        for i in self.list_of_peaks[gene_name][e]['peaks']:
                            plt.axvline(i-self.five_prime_flank, color='blue')
                    if self.print_valleys == True:
                        for i in self.list_of_peaks[gene_name][e]['valleys']:
                            plt.axvline(i-self.five_prime_flank, color='red')
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.2, color='green')

                # #mark Ts
                #     new_data = self.data[gene_name][self.data[gene_name].nucleotide == 'T']
                #     for n in list(new_data["position"]):
                #         plt.axvline(n, color='grey')

                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e][RT_begin:])
                    x_array = np.array(self.data[gene_name]['position'][RT_begin:])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array, color='grey')
                    ax2 = plt.twinx()
                    ax2.set_ylabel("no. of T per 9 nt")
                    plt.plot(self.data[gene_name]['position'], self.data[gene_name]['Trichness'], color="#800000")
                    plt.ylim(ymax=15)
                    # plt.text(210-self.five_prime_flank,max(self.data[gene_name][e])-150,'RT='+str(round(self.genes[gene_name]['RT'][e],3)))
                except KeyError:
                    plt.xlabel("NO READS")
                if plot_no == 6:
                    fig_no += 1
                    plt.savefig(self.prefix+'Trichness_'+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no)+'.png')
                    plt.clf()
                    plot_no = 0
            if plot_no > 0:
                plt.savefig(self.prefix+'Trichness_'+e+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'_fig_'+str(fig_no+1)+'.png')
                plt.clf()
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
            fig = plt.figure(figsize=(12, 9), facecolor='w', edgecolor='k')
            gene_name   =   self.id_to_names[i_gene_id]
            gene_length =   self.genes[gene_name]['gene_length']
            RT_begin    =   self.five_prime_flank + gene_length + self.readthrough_start
            plot_no = 0
            for e in self.experiments:
                plot_no += 1
                fig.add_subplot(subplot_layout[0],subplot_layout[1],plot_no)
                fig.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name]['position'], self.data[gene_name][e])
                    if self.print_peaks == True:
                        for i in self.list_of_peaks[gene_name][e]['peaks']:
                            plt.axvline(i-self.five_prime_flank, color='blue')
                    if self.print_valleys == True:
                        for i in self.list_of_peaks[gene_name][e]['valleys']:
                            plt.axvline(i-self.five_prime_flank, color='red')
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.2, color='green')
                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e][RT_begin:])
                    x_array = np.array(self.data[gene_name]['position'][RT_begin:])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                    # plt.text(0.5,0.5,'RT='+str(round(self.genes[gene_name]['RT'][e],3)))
                except KeyError:
                    plt.xlabel("NO READS")
                if plot_no == len(self.experiments):
                    plot_no = 0
                    plt.savefig(self.prefix+i_gene_id+'_l'+str(self.lookahead)+'_t'+str(self.hits_threshold)+'.png', dpi=150)
                    plt.clf()
        return True

# print plots, gene after gene,
    def fig_gene_after_gene(self):
        print '# Plotting gene after gene.'
        for e in self.experiments:
            fig = plt.figure(figsize=(12, 9), dpi=200, facecolor='w', edgecolor='k')
            fig_no = 0
            plot_no = 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                gene_length =   self.genes[gene_name]['gene_length']
                RT_begin    =   self.five_prime_flank + gene_length + self.readthrough_start
                plot_no += 1
                fig.add_subplot(3, 2, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name]['position'], self.data[gene_name][e])
                    if self.print_peaks == True:
                        for i in self.list_of_peaks[gene_name][e]['peaks']:
                            plt.axvline(i-self.five_prime_flank, color='blue')
                    if self.print_valleys == True:
                        for i in self.list_of_peaks[gene_name][e]['valleys']:
                            plt.axvline(i-self.five_prime_flank, color='red')
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.2, color='green')
                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e][RT_begin:])
                    x_array = np.array(self.data[gene_name]['position'][RT_begin:])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                    # plt.text(210-self.five_prime_flank,max(self.data[gene_name][e])-150,'RT='+str(round(self.genes[gene_name]['RT'][e],3)))
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

# print plots, gene after gene,
    def mark_T(self, anti_plot=False):
        print '# Plotting gene after gene marking each T.'
        for e in self.experiments:
            fig = plt.figure(figsize=(12, 9), dpi=200, facecolor='w', edgecolor='k')
            fig_no = 0
            plot_no = 0
            for i_gene_id in self.genes_id_list:
                gene_name = self.id_to_names[i_gene_id]
                gene_length =   self.genes[gene_name]['gene_length']
                RT_begin    =   self.five_prime_flank + gene_length + self.readthrough_start
                plot_no += 1
                fig.add_subplot(3, 2, plot_no)
                plt.tight_layout()
                plt.title(e)
                plt.ylabel("no. of reads")
                try:
                    plt.plot(self.data[gene_name]['position'], self.data[gene_name][e])
                    if self.print_peaks == True:
                        for i in self.list_of_peaks[gene_name][e]['peaks']:
                            plt.axvline(i-self.five_prime_flank, color='blue')
                    if self.print_valleys == True:
                        for i in self.list_of_peaks[gene_name][e]['valleys']:
                            plt.axvline(i-self.five_prime_flank, color='red')
                    for i in self.genes[gene_name]['exons']:
                        plt.axvspan(i[0], i[1], alpha=0.2, color='green')

                    #marking all T nucleotides
                    new_data = self.data[gene_name][self.data[gene_name].nucleotide == 'T']
                    for n in list(new_data["position"]):
                        plt.axvline(n, color='grey')
                    if anti_plot == True:
                        new_data = self.data[gene_name][self.data[gene_name].nucleotide == ('G' or 'C')]
                        for n in list(new_data["position"]):
                            plt.axvline(n, color='red', alpha=0.7)

                # some trick figured out by Hywel to plot without problems
                    y_array = np.array(self.data[gene_name][e][RT_begin:])
                    x_array = np.array(self.data[gene_name]['position'][RT_begin:])
                    y_array[0] = 0
                    y_array[len(y_array)-1] = 0
                    plt.fill(x_array, y_array)
                    plt.xlabel('ID: '+i_gene_id+', Name: '+gene_name)
                    # plt.text(210-self.five_prime_flank,max(self.data[gene_name][e])-150,'RT='+str(round(self.genes[gene_name]['RT'][e],3)))
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
            fig = plt.figure(figsize=(12, 9), dpi=200, facecolor='w', edgecolor='k')
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
                        plt.axvspan(i[0], i[1], alpha=0.2, color='green')
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
                        plt.axvspan(i[0], i[1], alpha=0.2, color='green')
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
                            r_delta_H += self.rna_dna[pair][0]
                            r_delta_S += self.rna_dna[pair][1]
                            d_delta_H += self.dna_dna[pair][0]
                            d_delta_S += self.dna_dna[pair][1]
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