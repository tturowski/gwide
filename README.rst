=====
gwide
=====

Instalation
===========
for users: 
pip install --user mercurial git+git://github.com/tturowski/gwide.git (or skip --user mercurial to install globally)
add to your .profile or .bash_profile or .bash_login
PATH="$HOME/.local/bin:$PATH"

for developers:
clone folder
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

gwideHittable (-h for help) - Three options:
  - Calculate "correlations": Pearson, Spearman or Kendall Tau
  - "count" hittables for further analysis. Ideal to work with multiple experiments
  - Plot "piechart"s for hittable classes (plots are not ideal...)

gwidePlot
  - plot genome wide plots, 5' and 3' end aligned or aligned to choosen aligner (-o aligner) or aligneg to 3' end of read-through
  - plot genome wide ratio between different experiment
  - possibility to filter genes using -f option
  - printing *.csv tables to generate heatmaps using other software (i.e. GENE-E)
  - calculate p-value for a non-canonical termination sites 
  - making GTF files (for whole transcripts or only for extensions

Old description
===========


pyPiecharts_from_hittables.py, compareHittables.py and countHittables.py working with pyReadCounter.py output files *hittable_reads.txt

compareHittables.py - to calculate Pearson, Spearman and Kendall corelation coefficients between hittables - as many as you want.
ruffusCalculatePearsonCoefficient.py - runs compareHittables.py but work directly on *.novo files

ruffusCreateConcat.py - creates *.concat file from *.novo files taking name from name.novo as experiment name
gwidePlot.py - generates genome-wide plots from *.concat files


Note
====

This project has been set up using PyScaffold 2.4.4. For details and usage
information on PyScaffold see http://pyscaffold.readthedocs.org/.
