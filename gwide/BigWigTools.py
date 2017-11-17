import sys, os, warnings
import random, pyBigWig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gwide.profileAnalyser as pa

def BigWigCount(datasets=dict(), genes=dict(), normalization="pM"):
    '''Calculates coverage for individual genes for given experiments

    datasets : dict()
      Dictionary with exp name as key and pyBigWig object

    genes : dict()
      Dictionary with pyBigWig.chroms() dictionary

    normalization : str()
      Normalization option" '"pM", "pMK", None'. Default: '"pM"'

    Returns
    -------
    pd.DataFrame()
    '''
    output_df = pd.DataFrame()
    for name, exp in datasets.iteritems():
        scaling_factor = datasets[name].header()['sumData'] / 1000000.0  # for pM
        exp_dict = dict()

        if normalization == "pM":
            for gene_name, length in genes.iteritems():  # normalizing to coverage per milion
                exp_dict[gene_name] = sum(datasets[name].values(gene_name, 0, length)) / scaling_factor

        elif normalization == "pMK":
            for gene_name, length in genes.iteritems():  # normalizing to coverage per milion per kb
                kb_scaling_factor = length / 1000.0 #gene length in kb
                exp_dict[gene_name] = (sum(
                    datasets[name].values(gene_name, 0, length)) / scaling_factor) / kb_scaling_factor

        elif normalization == None:
            for gene_name, length in genes.iteritems():
                exp_dict[gene_name] = sum(datasets[name].values(gene_name, 0, length))

        s1 = pd.Series(data=exp_dict, name=name)
        output_df = output_df.append(s1)

    return output_df.T

def BigWigPileup(datasets=dict(), gene_name=str(), transcript_length=int(), normalization='ntotal', window=1,
                 median_only=False):
    '''Returns pd.DataFrame() for multiple experiments with calculated median or all statistic

    datasets : dict()
      Dictionary with exp name as a key and pyBigWig object

    gene_name : str()

    normalization : str()
      Normalization option" '"pM", "ntotal", "raw"'. Default: '"ntotal"'

    median_only : boolean
      Default = '"False"' .If '"True"' retruns only pd.Series() for median

    window : int()
      Default = 1. Window size for rolling average to smooth data

    Returns
    -------
    pd.DataFrame() or pd.Series()
    '''
    data_df = pd.DataFrame()

    # looping through experiments and extracting puleups
    for name, exp in datasets.iteritems():
        normalizator = {'raw': 1,
                        'pM': (datasets[name].header()['sumData'] / 1000000),  # per Million
                        'ntotal': sum(datasets[name].values(gene_name, 0, transcript_length))}
        s1 = pd.Series(datasets[name].values(gene_name, 0, transcript_length)) / normalizator[normalization]
        data_df[name] = s1

    # smoothing
    if window > 1:
        data_df = data_df.rolling(window, win_type='blackman', center=True).mean()

    # calculating output
    if median_only == True:
        return data_df.median(axis=1)
    else:
        output_df = pd.DataFrame()
        for function in ['mean', 'median', 'min', 'max']:
            output_df[function] = getattr(data_df, function)(axis=1)  # calculates using pandas function listed in []
        if len(data_df.columns) > 2:  # calculating quartiles only in more than two experiments
            output_df['q1'], output_df['q3'] = data_df.quantile(q=0.25, axis=1), data_df.quantile(q=0.75, axis=1)
        return output_df

def mRNAalign(df_genes=pd.DataFrame(), datasets=dict(), normalization='ntotal', window=1):
    '''Returns DataFrame() for multiple experiments with calculated median and etc.

    df_genes : DataFrame()
      DataFrame() containing genes to align

    datasets : dict()
      Dictionary with exp name as a key and pyBigWig object

    normalization : str()
      Normalization option" '"pM", "ntotal", "raw"'. Default: '"ntotal"'

    window : int()
      Default = 1. Window size for rolling average to smooth data

    Returns
    -------
    dict() containing DataFrames for different aligners
    '''
    # constants to subtract
    for_5end = df_genes.max()['5UTR']
    for_3end = df_genes.max()['5UTR_and_CDS']
    for_pA = df_genes.max()['len']

    # values to be added to align
    df_genes['for_5end'] = for_5end - df_genes['5UTR']
    df_genes['for_3end'] = for_3end - df_genes['5UTR_and_CDS']
    df_genes['for_pA'] = for_pA - df_genes['len']

    # filling dataframes
    df_very5end = pd.DataFrame()
    df_5end = pd.DataFrame()
    df_3end = pd.DataFrame()
    df_pA = pd.DataFrame()
    for gene_name in df_genes.index.values:
        aligner = df_genes[df_genes.index == gene_name]  # details for the gene_name
        aligner = aligner.astype(int)
        s1_pileup = BigWigPileup(datasets, gene_name, int(aligner['len'] + 10.0), normalization, window,
                                 median_only=True)
        # very 5' end
        df_very5end[gene_name] = s1_pileup
        # 5' end
        to_add = pd.Series([0] * aligner['for_5end'][0])  # lenght to be added to align
        s2 = pd.concat([to_add, s1_pileup], ignore_index=True)
        df_5end[gene_name] = s2
        # 3' end
        to_add = pd.Series([0] * aligner['for_3end'][0])  # lenght to be added to align
        s2 = pd.concat([to_add, s1_pileup], ignore_index=True)
        df_3end[gene_name] = s2
        # pA site
        to_add = pd.Series([0] * aligner['for_pA'][0])  # lenght to be added to align
        s2 = pd.concat([to_add, s1_pileup], ignore_index=True)
        df_pA[gene_name] = s2

    return {"very 5' end": df_very5end, "5' end": df_5end, "3' end": df_3end, "pA site": df_pA}


def metaplot(df, df2=pd.DataFrame(), h_line=int(), title=str(), color='black', lc='red', figsize=(10, 6), dpi=200,
             center=True):
    '''Creates metaplot

    df : DataFrame()

    df2 : DataFrame()
        Dataframe with GC%

    h_line : int()
        Line visualazing aligning position

    title : str()

    color : str()

    lc    : str()
        horozontal line color
    figsize : touple()

    dpi  : int()

    center : boolean
        detault = True; center plot 1000 bp arount aligning site
    '''
    df = df.T.sum()  # plotting metaplot
    if center == True:  # options to center 1000 around aligning site
        if h_line - 1000 < 0:
            df = df[0:1000]
        elif h_line - 500 > 0 and h_line + 500 < len(df):
            df = df[int(h_line) - 500:int(h_line) + 500]
        else:
            df = df[-1000:]
    sns.set(font_scale=1.5)
    sns.set_style('white')
    sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads')
    #     ax1.set_ylim(ylim)
    ax1.plot(df, color=color)
    ax1.axvline(h_line, color=lc)

    if not df2.empty:
        df2 = df2.T.mean()
        if center == True:  # options to center 1000 around aligning site
            if h_line - 1000 < 0:
                df2 = df2[0:1000]
            elif h_line - 500 > 0 and h_line + 500 < len(df):
                df2 = df2[int(h_line) - 500:int(h_line) + 500]
            else:
                df2 = df2[-1000:]
        ax2 = ax1.twinx()
        ax2.plot(df2, color='orange')
        ax2.set_ylabel('GC%', color='orange')
        ax2.set_yticks(np.arange(0, 1, 0.1))
        ax2.set_ylim(0, 1)
        for tl in ax2.get_yticklabels():
            tl.set_color('orange')
    ax1.grid()
    ax1.legend()
    plt.show()


def featureBinCollect(data=pd.Series(), df_genes=pd.DataFrame(), gene_name=str(), bins=[1, 10, 1], warn=False):
    '''Calculates coverage districution for mRNA elements (UTRs and CDS)

    data : Series()

    df_genes : DataFrame()
      DataFrame() containig gene_details
    gene_name : str()

    bins : list()
        Number of bins for each feature [5UTR, CDS, 3UTR]. Default = '[1,10,1]'

    Returns
    -------
    Series()'''

    details = df_genes[df_genes.index == gene_name]  # details for the given gene
    if len(data) != int(details['len']):
        len_difference = len(data) - int(details['len'])
        data = data[:-len_difference]  # removing polyA when applicable
        if warn == True: warnings.warn("Data " + str(len_difference) + " nt longer than expected from given details.",
                                       UserWarning)
    normalized_data = data / data.sum()  # normalize data

    def binSum(data, bins):
        splited = np.array_split(np.array(data), bins)  # splitting into an even chunks
        return list(pd.DataFrame(splited).sum(1))  # returns sum

    output = list()
    # 5' UTR
    a = normalized_data[:int(details['5UTR'])]
    output += binSum(a, bins[0])
    # CDS
    b = normalized_data[int(details['5UTR']):int(details['5UTR']) + int(details['CDS'])]
    output += binSum(b, bins[1])
    # 3' UTR
    c = normalized_data[int(details['5UTR']) + int(details['CDS']):]
    output += binSum(c, bins[2])
    return pd.Series(output)


def codons(data=pd.Series(), df_genes=pd.DataFrame(), df_sequences=pd.DataFrame(), gene_name=str()):
    '''Annotates data with nucleotide, codone and aminoacid

    data : Series()

    df_genes : DataFrame()
      DataFrame() containig gene_details
    gene_name : str()

    df_sequences : DataFrame()
        sequences

    Returns
    -------
    DataFrame()'''

    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

    def triple_codons(sequence):
        output_list = list()
        for i in range(0, len(sequence), 3):
            output_list.append(sequence[i:i + 3])
            output_list.append(sequence[i:i + 3])
            output_list.append(sequence[i:i + 3])
        return output_list

    sequence = df_sequences[df_sequences.index == gene_name]['sequence'][0]  # sequence
    details = df_genes[df_genes.index == gene_name]  # details for the given gene

    output_df = pd.DataFrame()
    output_df['CDS'] = data[int(details['5UTR']):int(details['5UTR']) + int(details['CDS'])]
    sequence = sequence[int(details['5UTR']):int(details['5UTR']) + int(details['CDS'])]  # +3 includes stop codone
    output_df['sequence'] = list(sequence)
    if len(output_df) % 3 != 0: warnings.warn("Sequence length do not divided by 3", Warning)
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]  # codons generator
    output_df['codon'] = triple_codons(sequence)  # codons generator x3
    output_df['aminoacid'] = output_df['codon'].replace(codontable)  # codon to aminoacid
    return output_df