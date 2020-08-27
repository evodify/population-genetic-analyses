# Population Genetic Analyses

This repository contains a set of tools to perform various population genetic analyses on genomic data.

[ABBA-BABA](ABBA-BABA) folder contains scripts to run the ABBA-BABA test of admixture.

[calculate_AveragePerWindow.py](calculate_AveragePerWindow.py) calculates average values for non-overlapping sliding windows.

[calculate_WindowAverage_perSNP.py](calculate_WindowAverage_perSNP.py) calculates average values for a window centered on each SNP.

[calculate_MedianPerWindow.py](calculate_AveragePerWindow.py) calculates median values per sliding window.

[calculate_FixedHetero_PerWindow.py](calculate_FixedHetero_PerWindow.py) calculates fixed heterozygosity (exist in hybrids) with the sliding window approach.

[calculate_Hetero_PerWindow.py](calculate_Hetero_PerWindow.py) calculates heterozygosity with the sliding window approach.

[calculate_iHSproportion.py](calculate_iHSproportion.py) calculates fractions of SNPs with [iHS](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072) values above 2.0 over genomic windows of specified size.

[calculate_PositionMean.py](calculate_PositionMean.py) calculates mean per position of a table with numbers.

[calculate_Total-Hetero.py](calculate_Total-Hetero.py) calculates total heterozygosity.

[calculate_TotalMean.py](calculate_TotalMean.py) calculates total mean of a table with numbers.

[DAF.py](DAF.py) calculates change in derived allele frequency (deltaDAF).

[mutMatrix.py](mutMatrix.py) estimates a mutation matrix required for [LDhelmet](http://dx.doi.org/10.1371/journal.pgen.1003090).

[SFS.py](SFS.py) reconstructs the site frequency spectrum for a given set of samples.

**DISCLAIMER:** USE THESE SCRIPTS AT YOUR OWN RISK. I MAKE NO WARRANTIES THAT THESE SCRIPTS ARE BUG-FREE, COMPLETE, AND UP-TO-DATE. I AM NOT LIABLE FOR ANY LOSSES IN CONNECTION WITH THE USE OF THESE SCRIPTS.
