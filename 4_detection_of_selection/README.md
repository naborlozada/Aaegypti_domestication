## Selection

This directory contains different script files that are divided based on the type of analysis used to detect Selection (directory names): 1) PCadapt, 2) RAiSD, 3) MKT test and DoS statistics.

##
**PCadapt** 

Directory: PCadapt.

* **`outliers_bestK_pcadapt.R`:** R script to identify the best K used to identify outlier SNPs across populations using PCadapt.
* **`outliers_pcadapt.R`:** R script to identify the total set of outlier of SNPs of a population using PCadapt.


##
**RAiSD**

Directory: RAiSD.

These scripts are divided in those used to detect outliers based on (1) the mahalanobis distances and (2) top 0.05% and 1% tails of the `mu statistics` distribution from the RAiSD ouput.

* **`RAiSD_hard_sweeps.md`:** Descriptive steps to do RAiSD analysis at local or global level in our populations by chromosome.

(1) Outliers based on the mahalanobis distances, using as threshold FDR of 5% (alpha=0.05). 
* **`outliers_RAiSD.mahalanobis_dist.R`:** R script to identify the total set of outliers in each RAiSD analysis at local or global level.

(2) Outliers based on the top tails (0.05% and 1%) of the original raw `mu statistics` distribution.
* **`selective_sweeps.raisd.get_tail_tops_1_n_0.05.R` :** R script that performs two main cutoffs on the output and creates, for each cutoff, two output file containing all analyzed windowed sites and two BED outfiles to be used for mapping against the reference genome of *Ae. aegypti* (AaegL5).
* **`gene_distribution_across_Africa_and_Out_of_Africa.LOCAL_RAiSD_analysis.pl` :** PERL script taht takes the filtered output files of each population parsed from the above R script and analyze each gene case under selection whether it is present in each population. Next, it detects cases that are shared and/or population specific among all populations, but particularly focused on distinguishing cases in Out of Africa specific, but also shows other cases: shared between Africa and Out of Africa, and Africa population specific.
    
##
**Detection of Selection using the McDonald and Kreitman test (MKT) and the Direction of Selection (DoS)**

Directory: MKT_and_DoS.

Description of each steps and scripting codes names (bash, R, PERL) is in the page: `Dection_of_Selection_with_MKT.md`. 
