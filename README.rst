=====
gwide
=====



Description
===========

ploting of genome-wide plots pyCRAC downstream tool

gwideHittable (-h for help) - Three options:
  - Calculate "correlations": Pearson, Spearman or Kendall Tau
  - "count" hittables for further analysis. Ideal to work with multiple experiments
  - Plot "piechart"s for hittable classes (plots are not ideal...)

gwidePlot



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
