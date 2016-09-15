import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def compare1toRef(dataset=pd.Series(), ranges='mm', heatmap=False, relative=False,
                    reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    '''Takes series and compare this with reference DataFrame'''
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
        return differences_df['rel_diff']  # return Serie
    elif heatmap == True:
        return differences_df['diff']  # return Serie
    if heatmap == False:
        return_df['exp'] = dataset
        for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max']: return_df[f] = differences_df[f]
        return return_df  # return Dataframe

def plot_heatmap(df, title='Heat map of differences between dataset and reference plot for RDN37-1', vmin=None,
                 vmax=None):
    fig, ax = plt.subplots(figsize=(20, 10))
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
    working_df, result_df = pd.DataFrame(), pd.DataFrame()
    print "Experiments:"
    for f in [d for d in list(input_df.columns.values) if any(i in d for i in let_in) and any(o not in d for o in let_out)]:
        print f
        working_df[f]=input_df[f]
    for function in ['mean', 'median', 'min', 'max']: result_df[function]=getattr(working_df, function)(axis=1)
    if len(working_df.columns) > 2: #calculating quartiles only in more than two experiments
        result_df['q1'], result_df['q3'] = working_df.quantile(q=0.25, axis=1),working_df.quantile(q=0.75, axis=1),
    return result_df

def plot_to_compare(dataset=pd.DataFrame(), label=str(), start=None, stop=None, figsize=(15,8), lim=0.01 ,reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    reference = pd.read_csv(reference, index_col=0) # reading reference
    dataset, s2 = dataset[start:stop], reference[start:stop] # prepating datasets
    #plotting
    fig, ax1 = plt.subplots(figsize=figsize)
    plt.axhline(0, color='red')
    if len(dataset.columns) == 4: #if only two experiments
        ax1.plot(dataset.index, dataset['mean'], 'black', label=label)
        ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.3)
    else: #if more than two experiments
        ax1.plot(dataset.index, dataset['median'], 'black', label=label)
        ax1.fill_between(dataset.index, dataset['q1'], dataset['q3'], color='black', alpha=0.2)
        ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.07)
    ax1.set_xlabel('position')
    ax1.set_ylim(0,lim)
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('fraction of reads '+label, color='black')
    for tl in ax1.get_yticklabels():
        tl.set_color('black')

    ax1.plot(s2.index, s2['median'], 'green', label='reference RDN37-1')
    ax1.fill_between(s2.index, s2['q1'], s2['q3'], color='green', alpha=0.2)
    ax1.fill_between(s2.index, s2['min'], s2['max'], color='green', alpha=0.07)
    ax1.legend()

def compareMoretoRef(dataset=pd.DataFrame(), ranges='mm',
                     reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    '''Takes series and compare this with reference DataFrame'''
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
              plot_ranges=True, figsize=(15, 8), lim=0.01,
              reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
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
            ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.3)
        else:  # if more than two experiments
            ax1.fill_between(dataset.index, dataset['q1'], dataset['q3'], color='black', alpha=0.2)
            ax1.fill_between(dataset.index, dataset['min'], dataset['max'], color='black', alpha=0.07)
        ax1.fill_between(s2.index, s2['q1'], s2['q3'], color='green', alpha=0.2)
        ax1.fill_between(s2.index, s2['min'], s2['max'], color='green', alpha=0.07)
    ax1.set_ylim(0, lim)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads ' + label, color='black')
    plt.legend()
