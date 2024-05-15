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
**pN/pS ratio**

(1) Based on only segregating sites of Ae. aegypti's populations:
* **`pnps_ratio.R`:** R script to calculate the pn/pn ratio of SNPs of a population.

(2) Based on the divergence between two sister species: Ae. aegypti and Ae. albopictus using YN model from PAML. 
* **`detection_of_selection_PAML.md`** Summary protocol steps and the script used to calculate the pN/pS ratio for each protein coding gene (~500,000) in each single population.
* **`merge_orthologues_genes.pl`** PERL script ot join two orthologs genes in a single fasta file.
* **`GENEID_yn00.ctl`** Control infile from PAML YN model (YN_model.ctl). 
* **`get_pnps_ratio_per_poplation.pl`** PERL script that extract the pN/pS statistics from PAML analysis and sort them for each single gene present in a population, as well as defined the type of selection based on the threshold defined in `Methods section` of this paper.


##
Comment on the **Update for the rebuttal 2024**\
Two additional scripts that were used to performed the outlier detection based on the top 0.05% and 1% tail distribution of the RAiSD results were added in this directory:\
\
**RAiSD**
* R script: `selective_sweeps.raisd.LOCAL.get_tail_tops_1_n_0.05.R`.
* PERL script: `gene_distribution_across_Africa_and_Out_of_Africa.LOCAL_RAiSD_analysis.pl`.

**pN/pS ratio**
* `detection_of_selection_PAML.md`
* `merge_orthologues_genes.pl`
* `GENEID_yn00.ctl`
* `get_pnps_ratio_per_poplation.pl`

