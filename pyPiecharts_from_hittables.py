#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
__copyright__	= "Copyright 2015"
__version__		= "0.2"
__credits__		= ["Tomasz Turowski"]
__email__		= "twturowski@gmail.com"
__status__		= "Production"

import os, argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Usage: Plots piecharts for all *hittable_reads.txt in current directory')
parser.add_argument("-s", "--single", dest="print_single", help="Print hittables in single files",
                     action="store_true", default=False)
args = parser.parse_args()

experiments = ["_".join(f.split('_')[0:3]) for f in os.listdir(os.getcwd()) if f.endswith("hittable_reads.txt")]

#initiating DataFrame
data = pd.DataFrame(columns=[['group']+['legend']+experiments]) # initialize Pandas DataFrame
data = data.set_index(['group'])
general = dict()

#filling DataFrame
for hittable in os.listdir(os.getcwd()):
    if hittable.endswith("hittable_reads.txt"):
        print 'Reading file: '+hittable
        name = "_".join(hittable.split('_')[0:3]) #take fist three elements of file name as experiment name
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

if args.print_single == False:
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
        # plt.legend(data.index, loc=4)
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

elif args.print_single == True:
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