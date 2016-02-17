=====
gwide
=====
Python package for downstream analysis of CRAC data.

Instalation
===========
for users: 
pip install --user mercurial git+git://github.com/tturowski/gwide.git (or skip --user mercurial to install globally);
add to your .profile or .bash_profile or .bash_login:
PATH="$HOME/.local/bin:$PATH"

for developers:
clone folder and;
pip install --editable .

Setup
===========
It's possible to create $HOME/bin/default.aml file with paths to follwing files:
#default paths to genome files
GTF_PATH : *.gtf
FASTA_PATH : *.fasta
TAB_PATH : *.tab

Description
===========

ploting of genome-wide plots pyCRAC downstream tool

gwideHittable (-h for help, input *.hittable) - Three options:
  - Calculate "correlations": Pearson, Spearman or Kendall Tau
  - "count" hittables for further analysis. Ideal to work with multiple experiments
  - Plot "piechart"s for hittable classes (plots are not ideal...)

novo2concat.py (-h for help):
  - copy *.novo files to new folder
  - prefix from *.novo file is used as experiment name

gwidePlot (-h for help, input *.concat):
  - plot genome wide plots, 5' and 3' end aligned or aligned to choosen aligner (-o aligner) or aligneg to 3' end of read-through
  - plot genome wide ratio between different experiment
  - possibility to filter genes using -f option
  - printing *.csv tables to generate heatmaps using other software (i.e. GENE-E)
  - calculate p-value for a non-canonical termination sites 
  - making GTF files (for whole transcripts or only for extensions

gwidetRNA (-h for help, input *.concat):
  - plot single tRNA plots: multiple experiments per page, one experiment per page, one experiment under another, mark A and B boxes
  - plot ratio between experiments
  - plot under nucleotides resolution
  - calculate dG RNA_DNA/DNA_DNA for each valley or last 20 nt after tRNA gene
  - save tab-deliminated file with each plot calculations

Old description
===========
ruffusCalculatePearsonCoefficient.py - runs compareHittables.py but work directly on *.novo files

Note
====
This project has been set up using PyScaffold 2.4.4. For details and usage
information on PyScaffold see http://pyscaffold.readthedocs.org/.
