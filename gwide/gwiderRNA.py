#!/usr/bin/env python
from optparse import OptionParser
import select
from gwide.Classes.rRNAFromConcat import *
import gwide.methods as gtm


""" Script working with concat file generated by pileupsToConcat.py script. Can work on stdin from that script. Read concat file and according to options.
    Print rRNA genes in order of experiments (3 per page)"""

def rRNA():
    usage = "Usage: create pileups with pyPileup (pyCRAC package) then in directory containing pileup files type run i.e.:"+"\n"+ \
            "cat file.concat | gwiderRNA.py or gwiderRNA.py -i file.concat"
    parser = OptionParser(usage=usage)
    parser.add_option("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     metavar="FILE", default=None)
    parser.add_option("-i", dest="input_file", help="Provide the path to your concat file.",
                     metavar="FILE", default=None)
    parser.add_option("--5flank", dest="five_prime_flank", type="int", help="Set up 5 prime flank in pileup file. Default = 1000", default=1000)
    parser.add_option("--3flank", dest="three_prime_flank", type="int", help="Set up 3 prime flank in pileup file. Default = 1000", default=1000)
    parser.add_option("-l", "--lookahead", dest="lookahead", type="int", help="Set up lookahead parameter for pypeaks function. Default = 20", default=20)
    parser.add_option("-t", "--hits_threshold", dest="hits_threshold", type="int", help="Set up threshold for pileup. Default 100 reads",
                      default=100)
    # parser.add_option("-r", "--readthrough", dest="readthrough", type="int", help="Set up when readthrough should start countin. Default: 0",
    #                   default=0)
    parser.add_option("-p", "--prefix", dest="out_prefix", type="str", help="Prefix for output files. Default to standard output", default=None)
    parser.add_option("--peaks", dest="print_peaks", action="store_true", help="Add into command line if you want to print peaks on plots. Default: False",
                      default=False)
    parser.add_option("-o", "--output", dest="output_files", choices=["std", "ratio", "single", "correlations", "ratio_smooth"], help="Select from following options: (1) std - RDN37-1; experiment after experimen ;"+'\n'
                                                                                                   "(2)ratio - ratio for -a divided by -b; (3)single - plot RDN37-1 plots 1 per page; (4) correlations - calculate correlations for different experiments; (5)ratio_smooth - ratio for -a divided by -b", default="std")
    parser.add_option("-a", dest="to_divide", type="str", help="experiment to divide by -b", default=None)
    parser.add_option("-b", dest="divisor", type="str", help="experiment being divisor for -a", default=None)
    parser.add_option("-n", "--normalized", dest="normalized", action="store_true", help="Use when you want to work on data normalized reads per Milion? Default: False", default=False)
    (options, args) = parser.parse_args()

    gtf_file = gtm.getGTF(options.gtf_file)

    if options.out_prefix:
        prefix = options.out_prefix+'_'
    else:
        prefix = str()

    if options.output_files == "ratio":
        options.normalized = True

    data = rRNAFromConcat(gtf_file=gtf_file, five_prime_flank=options.five_prime_flank, print_peaks=options.print_peaks,
                             three_prime_flank=options.three_prime_flank, hits_threshold=options.hits_threshold, lookahead=options.lookahead, prefix=prefix, normalized=options.normalized)
    data.read_csv(options.input_file)
    data.slice_data()
    if options.print_peaks == True:
        data.find_peaks()
    if options.output_files == "std":
        data.print_rRNA()   # RDN37 should be prepared with 1000 nt flanks
    if options.output_files == "single":
        data.single_rRNA()   # RDN37 should be prepared with 1000 nt flanks
    if options.output_files == "ratio":
        # data.fig_ratio(options.to_divide, options.divisor)  # plots ratio to_divide/divisor
        data.fig_log2ratio(options.to_divide, options.divisor)  # plots log2 ratio to_divide/divisor
    if options.output_files == "ratio_smooth":
        data.fig_smoothlog2ratio(options.to_divide, options.divisor)  # plots log2 ratio to_divide/divisor using smoothed data
    if options.output_files == "correlations":
        data.correlations()
    print '# Done.'