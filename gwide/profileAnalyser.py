import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_from_csv(csv_path=str(),title=None, start=None, stop=None,figsize=(15,7),ylim=(None,0.01), color='green', h_lines=list(), lc='red'):
    '''Function creates plot similar to box plot: median, 2 and 3 quartiles and min-max range'''
    reference = pd.read_csv(csv_path, index_col=0).drop('nucleotide', 1)
    s2 = reference[start:stop]
    #plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads')
    ax1.set_ylim(ylim)
    ax1.plot(s2.index, s2['median'], color=color)
    if set(['q1','q3']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index, s2['q1'], s2['q3'], label='range (q2-q3)', color=color, alpha=0.2)
    if set(['min','max']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index, s2['min'], s2['max'], label='range (min-max)', color=color, alpha=0.07)
    for i in [i for i in h_lines if i in range(start,stop)]: ax1.axvline(i, color=lc)
    ax1.legend()

def plot_to_compare(dataset=pd.DataFrame(), label=str(), start=None, stop=None, figsize=(15,8), ylim=(None,0.01), h_lines=list() ,reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    '''Function creates two plots similar to box plot: median, 2 and 3 quartiles and min-max range
        ax1 - dataset given by user
        ax2 - reference dataset'''
    reference = pd.read_csv(reference, index_col=0) # reading reference
    dataset, s2 = dataset[start:stop], reference[start:stop] # prepating datasets
    #plotting
    fig, ax1 = plt.subplots(figsize=figsize)
    plt.axhline(0, color='red')
    if len(dataset.columns) == 4: #if only two experiments
        ax1.plot(dataset.index, dataset['mean'], 'black', label=label)
        ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.3, label='range (min-max)')
    else: #if more than two experiments
        ax1.plot(dataset.index, dataset['median'], 'black', label=label)
        ax1.fill_between(dataset.index, dataset['q1'], dataset['q3'], color='black', alpha=0.2, label='range (q2-q3)')
        ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.07, label='range (min-max)')
    ax1.set_xlabel('position')
    ax1.set_ylim(ylim)
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('fraction of reads '+label, color='black')
    for tl in ax1.get_yticklabels():
        tl.set_color('black')

    ax1.plot(s2.index, s2['median'], 'green', label='reference RDN37-1')
    ax1.fill_between(s2.index, s2['q1'], s2['q3'], color='green', alpha=0.2, label='range (q2-q3)')
    ax1.fill_between(s2.index, s2['min'], s2['max'], color='green', alpha=0.07, label='range (min-max)')
    for i in [i for i in h_lines if i in range(start, stop)]: ax1.axvline(i, color='red')
    ax1.legend()

def compare1toRef(dataset=pd.Series(), ranges='mm', heatmap=False, relative=False,
                    reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    '''Takes series and compare this with reference DataFrame, as a result gives:
    (a) if heatmap = False: Dataframe with(reference_above_experiment minimum etc.): rae_min, rae_max, ear_min, ear_max
    (b) if heatmap = True: Series of differences to plot heatmap
    relative : only for heatmap, recalculates differences according to the peak size. Warning: negative values are in range -1 to 0
    but positive are from 0 to values higher than 1 '''
    ranges_dict = {'mm': ['min', 'max'], 'qq': ['q1', 'q3']}

    # preparing dataframe and reference
    differences_df, return_df = pd.DataFrame(), pd.DataFrame()
    reference = pd.read_csv(reference, index_col=0)
    differences_df['exp'] = dataset
    differences_df['ref_min'] = reference[ranges_dict[ranges][0]]  # choosing q1 or min
    differences_df['ref_max'] = reference[ranges_dict[ranges][1]]  # choosing q3 or max

    ## finding differences (only positive value indicate difference)
    differences_df['ref_above_exp'] = differences_df['ref_min'] - differences_df['exp']
    differences_df['exp_above_ref'] = differences_df['exp'] - differences_df['ref_max']

    ## filling differences
    for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max', 'diff_ear', 'diff_rae', 'diff']: differences_df[f] = 0
    differences_df['rae_min'][differences_df['ref_above_exp'] > 0] = differences_df['exp']
    differences_df['rae_max'][differences_df['ref_above_exp'] > 0] = differences_df['ref_min']
    differences_df['ear_min'][differences_df['exp_above_ref'] > 0] = differences_df['ref_max']
    differences_df['ear_max'][differences_df['exp_above_ref'] > 0] = differences_df['exp']

    # combining differences for heatmap
    differences_df['diff_ear'] = differences_df['ear_max'] - differences_df['ear_min']
    differences_df['diff_rae'] = differences_df['rae_min'] - differences_df['rae_max']
    differences_df['diff'][differences_df['diff_ear'] > 0] = differences_df['diff_ear']
    differences_df['diff'][differences_df['diff_rae'] < 0] = differences_df['diff_rae']

    if heatmap == True and relative == True:
        differences_df['ref_median'] = reference['median']
        differences_df['rel_diff'] = differences_df['diff'] / differences_df['ref_median']
        return differences_df['rel_diff']  # return Series
    elif heatmap == True:
        return differences_df['diff']  # return Series
    if heatmap == False:
        return_df['exp'] = dataset
        for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max']: return_df[f] = differences_df[f]
        return return_df  # return Dataframe

def plot_heatmap(df=pd.DataFrame(), title='Heat map of differences between dataset and reference plot for RDN37-1', vmin=None,
                 vmax=None, figsize=(20,10)):
    '''plot heat map of differences, from dataframe generated by compare1toRef(dataset, heatmap=True) function'''
    fig, ax = plt.subplots(figsize=figsize)
    if not vmin:
        vmin = -np.absolute(df.max().median())
    if not vmax:
        vmax = df.max().median()
    heatmap = ax.pcolormesh(df.transpose(), cmap='seismic', vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(list(df.columns.values))) + 0.5, minor=False)
    ax.set_yticklabels(list(df.columns.values), minor=False)
    fig.colorbar(heatmap)
    ax.set_title(title)

def filter_df(input_df=pd.DataFrame(), let_in=[''], let_out=['wont_find_this_string']):
    '''using input dataframe, calculates mean, median etc. for all experiments "let in"'''
    working_df, result_df = pd.DataFrame(), pd.DataFrame()
    print "Experiments:"
    for f in [d for d in list(input_df.columns.values) if any(i in d for i in let_in) and any(o not in d for o in let_out)]:
        print f
        working_df[f]=input_df[f]
    for function in ['mean', 'median', 'min', 'max']: result_df[function]=getattr(working_df, function)(axis=1) #calculates using pandas function listed in []
    if len(working_df.columns) > 2: #calculating quartiles only in more than two experiments
        result_df['q1'], result_df['q3'] = working_df.quantile(q=0.25, axis=1),working_df.quantile(q=0.75, axis=1),
    return result_df

def compareMoretoRef(dataset=pd.DataFrame(), ranges='mm',
                     reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    '''Takes Dataframe created by filter_df and compare this with reference DataFrame'''
    ranges_dict = {'mm': ['min', 'max'], 'qq': ['q1', 'q3']}
    # preparing dataframe and reference
    differences_df, return_df = pd.DataFrame(), pd.DataFrame()
    reference = pd.read_csv(reference, index_col=0)
    if len(dataset.columns) == 4:  # if only two experiments
        differences_df['exp_min'] = dataset['min']
        differences_df['exp_max'] = dataset['max']
    else:  # if more than two exp
        differences_df['exp_min'] = dataset[ranges_dict[ranges][0]]  # choosing q1 or min
        differences_df['exp_max'] = dataset[ranges_dict[ranges][1]]  # choosing q3 or max
    differences_df['ref_min'] = reference[ranges_dict[ranges][0]]  # choosing q1 or min
    differences_df['ref_max'] = reference[ranges_dict[ranges][1]]  # choosing q3 or max

    ## finding differences (only positive value indicate difference)
    differences_df['ref_above_exp'] = differences_df['ref_min'] - differences_df['exp_max']
    differences_df['exp_above_ref'] = differences_df['exp_min'] - differences_df['ref_max']

    ## filling differences
    for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max']: differences_df[f] = 0
    differences_df['rae_min'][differences_df['ref_above_exp'] > 0] = differences_df['exp_max']
    differences_df['rae_max'][differences_df['ref_above_exp'] > 0] = differences_df['ref_min']
    differences_df['ear_min'][differences_df['exp_above_ref'] > 0] = differences_df['ref_max']
    differences_df['ear_max'][differences_df['exp_above_ref'] > 0] = differences_df['exp_min']
    #     return_df = dataset
    for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max']: return_df[f] = differences_df[f]
    return return_df  # returns Dataframe

def plot_diff(dataset=pd.DataFrame(), ranges='mm', label=str(), start=None, stop=None, plot_medians=True,
              plot_ranges=True, figsize=(15, 8), lim=0.01, h_lines=list(),
              reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    '''Plots differences calculated using compare1toRef or compareMoretoRef'''
    ranges_dict = {'mm': 'min-max', 'qq': 'q1-q3'}
    reference = pd.read_csv(reference, index_col=0)  # reading reference
    differences_df = compareMoretoRef(dataset=dataset, ranges=ranges)[start:stop]
    dataset, s2 = dataset[start:stop], reference[start:stop]  # prepating datasets
    # plotting
    fig, ax1 = plt.subplots(figsize=figsize)
    ax1.fill_between(differences_df.index, differences_df['ear_min'], differences_df['ear_max'], color='red',
                     where=(differences_df['ear_max'] > 0), label='increased pausing (' + ranges_dict[ranges] + ')')
    ax1.fill_between(differences_df.index, differences_df['rae_min'], differences_df['rae_max'], color='blue',
                     where=(differences_df['rae_max'] > 0), label='decreased pausing (' + ranges_dict[ranges] + ')')
    if plot_medians == True:
        ax1.plot(dataset.index, dataset['median'], 'black', label=label)
        ax1.plot(s2.index, s2['median'], 'green', label='reference RDN37-1')
    if plot_ranges == True:
        if len(dataset.columns) == 4:  # if only two experiments
            ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.3, label='min-max')
        else:  # if more than two experiments
            ax1.fill_between(dataset.index, dataset['q1'], dataset['q3'], color='black', alpha=0.2, label='q1-q3')
            ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.07, label='min=max')
        ax1.fill_between(s2.index, s2['q1'], s2['q3'], color='green', alpha=0.2, label='q1-q3')
        ax1.fill_between(s2.index, s2['min'], s2['max'], color='green', alpha=0.07, label='min=max')
    ax1.set_ylim(0, lim)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads ' + label, color='black')
    for i in [i for i in h_lines if i in range(start, stop)]: ax1.axvline(i, color='red')
    plt.legend()