## Selection

This directory contains different script files that are divided based on the type of analysis used to detect Selection: PCadapt, RAiSD, pN/pS ratio.

##
**PCadapt** 
* **`outliers_bestK_pcadapt.R`:** R script to identify the best K used to identify outlier SNPs across populations using PCadapt.
* **`outliers_pcadapt.R`:** R script to identify the total set of outlier of SNPs of a population using PCadapt.

##
**RAiSD**\
These scripts are divided in those used to detect outliers based on (1) the mahalanobis distances and (2) top 0.05% and 1% tails of the `mu statistics` distribution from the RAiSD ouput.

(1) Outliers based on the mahalanobis distances, using as threshold FDR of 5% (alpha=0.05). 
* **`RAiSD_hard_sweeps.md`:** Descriptive steps to do RAiSD analysis at local or global level in our samples.
* **`outliers_RAiSD.R`:** R script to identify the total set of outliers in each RAiSD analysis at local or global level.

(2) Outliers based on the top tails (0.05% and 1%) of the original raw `mu statistics` distribution.
* **`selective_sweeps.raisd.LOCAL.get_tail_tops_1_n_0.05.R` :** R script that performs two main cutoffs on the output and creates, for each cutoff, two output file containing all analyzed windowed sites and two BED outfiles to be used for mapping against the reference genome of *Ae. aegypti* (AaegL5).
* **`gene_distribution_across_Africa_and_Out_of_Africa.LOCAL_RAiSD_analysis.pl` :** PERL script taht takes the filtered output files of each population parsed from the above R script and analyze each gene case under selection whether it is present in each population. Next, it detects cases that are shared and/or population specific among all populations, but particularly focused on distinguishing cases in Out of Africa specific, but also shows other cases: shared between Africa and Out of Africa, and Africa population specific.
    
##
**Detection of Selection using the McDonald and Kreitman test (MKT) and the Direction of Selection (DoS)**

Summary steps description and code is in the page: `Dection_of_Selection_with_MKT.md`. 

##
Comment on the **Update for the rebuttal 2024**\
Two additional scripts that were used to performed the outlier detection based on the top 0.05% and 1% tail distribution of the RAiSD results were added in this directory:\
\
**RAiSD**
* R script: `selective_sweeps.raisd.LOCAL.get_tail_tops_1_n_0.05.R`.
* PERL script: `gene_distribution_across_Africa_and_Out_of_Africa.LOCAL_RAiSD_analysis.pl`.

**McDonald and Kreitman test (MKT) and the Direction of Selection (DoS)**
* `Dection_of_Selection_with_MKT.md`: contains all steps and a 'summary' code (as loops) used in a single protein coding genes poopulation, since all jobs were ran in parallel to optimize time.
