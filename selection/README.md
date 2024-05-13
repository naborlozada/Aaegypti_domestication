## Selection

This directory contains two files:
* **`outliers_bestK_pcadapt.R`:** R script to identify the best K used to identify outlier SNPs across populations using PCadapt.
* **`outliers_pcadapt.R`:** R script to identify the total set of outlier of SNPs of a population using PCadapt.
* **`pnps_ratio.R`:** R script to calculate the pn/pn ratio of SNPs of a population.
* **`RAiSD_hard_sweeps.md`:** Descriptive steps to do RAiSD analysis at local or global level in our samples.
* **`outliers_RAiSD.R`:** R script to identify the total set of outliers in each RAiSD analysis at local or global level.

Two additional scripts that were used to performed the outlier detection based on the top 0.05% and 1% tail distribution of the RAiSD results were added in this directory:
* R script (`selective_sweeps.raisd.LOCAL.get_tail_tops_1_n_0.05.R`): It performs two main cutoffs on the output and creates, for each cutoff, two output file containing all analyzed windowed sites and two BED outfiles to be used for mapping against the reference genome of *Ae. aegypti* (AaegL5).
* PERL script (`gene_distribution_across_Africa_and_Out_of_Africa.LOCAL_RAiSD_analysis.pl`). It takes the filtered output files of each population parsed from the above R script and analyze each gene case under selection whether it is present in each population. Next, it detects cases that are shared and/or population specific among all populations, but particularly focused on distinguishing cases in Out of Africa specific, but also shows other cases: shared between Africa and Out of Africa, and Africa population specific.
    
