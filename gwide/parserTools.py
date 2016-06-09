#!/usr/bin/env python
__author__ = 'Tomasz Turowski'
from pyCRAC.Parsers import GTF2
import sys, re
from optparse import *
from signal import signal, SIGPIPE, SIG_DFL
import gwide.methods as gtm

def getFastaSeqs():
    parser = OptionParser(usage="List of genes as std input and parameters")
    parser.add_option("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     type="str", default=None)
    parser.add_option("-f", "--fasta_file", dest="fasta_file", help="Provide the path to your fasta file.",
                     type="str", default=None)
    parser.add_option("-t", "--tab_file", dest="tab_file", help="Provide the path to your genom tab file.",
                     type="str", default=None)
    parser.add_option("-r", "--ranges", dest="ranges",
                     help="Provide ranges(flanks) for genes.",
                     type="int", default=0)
    parser.add_option("-a", "--5end", dest="five_end",
                     help="Set up 5` flank. If minus then print only 3` end. Python slicing [a:b] i.e. [200:401] - from 200 to 400; [-200:] - last 200; "
                          "[:-200] from begining till -200 before end",
                     type="int", default=None)
    parser.add_option("-b", "--3end", dest="three_end",
                     help="Set up 5` flank. If minus then print only 5` end. Python slicing [a:b]",
                     type="int", default=None)
    (options, args) = parser.parse_args()

    signal(SIGPIPE,SIG_DFL) # to manage with stdin and stdout
    #crating gtf object
    gtf = GTF2.Parse_GTF()
    gtf.read_GTF(gtm.getGTF(options.gtf_file))
    gtf.read_FASTA(gtm.getFASTA(options.fasta_file))
    gtf.read_TAB(gtm.getTAB(options.tab_file))

    for i in sys.stdin:
        gene_name = str(i.strip())
        genomic_seq = gtf.genomicSequence(gene_name, ranges=options.ranges)
        print '>'+gene_name
        print genomic_seq[options.five_end:options.three_end]+'\n'

def getGeneLength():
    parser = OptionParser(usage="usage: List of genes as std input")
    parser.add_option("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     type="str", default=None)
    (options, args) = parser.parse_args()

    signal(SIGPIPE,SIG_DFL)
    gtf = GTF2.Parse_GTF()
    gtf.read_GTF(gtm.getGTF(options.gtf_file))

    for i in sys.stdin:
        gene_name = str(i.strip())
        gene_length = gtf.geneLength(gene_name)
        print gene_name+"\t"+str(gene_length)


def getIdFromName():
    parser = OptionParser(usage="usage: List of genes as std input")
    parser.add_option("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     type="str", default=None)
    (options, args) = parser.parse_args()

    signal(SIGPIPE,SIG_DFL)

    gtf = GTF2.Parse_GTF()
    gtf.read_GTF(gtm.getGTF(options.gtf_file))

    for i in sys.stdin:
        gene_name = str(i.strip())
        gene_id = gtf.genes[gene_name]['gene_id']
        print gene_name+'\t'+gene_id

def getNameFromId():
    parser = OptionParser(usage="usage: List of genes as std input")
    parser.add_option("-g", "--gtf_file", dest="gtf_file", help="Provide the path to your gtf file.",
                     type="str", default=None)
    (options, args) = parser.parse_args()
    signal(SIGPIPE,SIG_DFL)

    gtf = GTF2.Parse_GTF()
    gtf.read_GTF(gtm.getGTF(options.gtf_file))

    id_to_gene = dict()
    for gene_name in gtf.genes:
        gene_id = gtf.genes[gene_name]['gene_id']
        id_to_gene[gene_id] = gene_name

    for i in sys.stdin:
        gene_id = str(i.strip())
        gene_name = id_to_gene[gene_id]
        print gene_name

def getGeneNamesFromGTF():
    parser = OptionParser(usage="getGenesNames; type usage: %prog [options] -f filename")
    files = OptionGroup(parser, "File input options")
    files.add_option("-f", "--input_file", dest="gtf_file",
                     help="Provide the path to your gtf data file. Default is standard input.",
                     type="str", default=None)
    files.add_option("-g", "--genes", dest="genes",
                     help="Which biotype of features to get: mRNA, tRNA, rRNA, snRNA, snoRNA",
                     type="str", default='tRNA')
    files.add_option("-i", "--introns", dest="introns",
                     help="Introns? both - not discriminate; int_cont -only intron containing; int_less - only int less",
                     choices=["both", "int_cont", "int_less"], default="both"),
    files.add_option("-o", "--output_file", dest="output_file",
                     help="Use this flag to provide an output file name. Default is standard output.", default=None)
    parser.add_option_group(files)
    (options, args) = parser.parse_args()

    ### By default, input and output are expected from the standard input or standard output.
    signal(SIGPIPE,SIG_DFL)
    outfile = sys.stdout
    if options.output_file:
        outfile = open(options.output_file, "w")

    gtf = GTF2.Parse_GTF()
    gtf.read_GTF(gtm.getGTF(options.gtf_file))
    names_list = list()

### for loop extracting tRNA names
    for line in open(gtm.getGTF(options.gtf_file), "r"):
        if not line.startswith('#'):
            line_elements = line.strip().split('\t')
            # assert len(line_elements) == 10, 'Unexpected number of elements found in gtf line: ' + line
            if str(line_elements[1]) == options.genes:
                try:
                    name = re.search("gene_name\s\"(.*?)\"", str(line_elements[8])).group(1)
                except:
                    pass
                #     name = re.search("gene_id\s\"(.*?)\"", str(line_elements[8])).group(1)
                if options.introns == "both":
                    if name not in names_list:
                        names_list.append(name)
                    # outfile.write(str(name) + '\n')
                elif options.introns == "int_cont":
                    if gtf.intronCoordinates(name):
                        if name not in names_list:
                            names_list.append(name)
                        # outfile.write(str(name) + '\n')
                elif options.introns == "int_less":
                    if not gtf.intronCoordinates(name):
                        if name not in names_list:
                            names_list.append(name)
                        # outfile.write(str(name) + '\n')

    outfile.write('\n'.join(names_list)+'\n')
    outfile.close()