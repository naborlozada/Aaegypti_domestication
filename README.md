# Aaegypti_domestication
---
Title:
**Adaptive genomic signatures of globally invasive popultions of the yellow fever mosquito *Aedes aegypti****

**Authors:** Alejandro Nabor Lozada-Chávez, Irma Lozada-Chávez, Niccolò Alfano, Umberto Palatini, Davide Sogliani, Samia Elfekih, Teshome Degefa, Maria V. Sharakhova, Athanase Badolo, Sriwichai Patchara, Mauricio Casas-Martinez, Bianca C. Carlos, Rebeca Carballar-Lejarazú, Luis Lambrechts, Jayme A. Souza-Neto & Mariangela Bonizzoni

**ACCEPTED MANUSCRIPT**

--- 

<ins>This is a final version of the repository with all main scripts used for the study and Supplementary Data files cited in the manuscript. If you use one or more of these scripts or data in a complete or partial form for your research, as well as any other data contained here, **please cite us**: </ins>

This repository **`Github` and `Zenodo`**:\
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10721413.svg)](https://doi.org/10.5281/zenodo.7863455)
   
   * Lozada-Chávez, A.N., Lozada-Chávez, I., Alfano, N., Palatini, U., Sogliani, D., Elfekih, D., Degefa, T., Sharakhova, M.V., Badolo, A., Patchara, S., Casas-Martinez, M., Carlos B.C., Carballar-Lejarazú, R., Lambrechts, L., Souza-Neto, J.A., & Bonizzoni M. Adaptive genomic signatures of globally invasive populations of the yellow fever mosquito Aedes aegypti. Zenodo. https://doi.org/10.5281/zenodo.7863455 (2024).

The preprint **bioRxiv** paper:
   * A. N. Lozada-Chávez et. al. 2023. [Molecular signature of domestication in the arboviral vector *Aedes aegypti*](https://doi.org/10.1101/2023.03.13.532092). bioRxiv. DOI: https://doi.org/10.1101/2023.03.13.532092 [submitted]. 


## Description

This repository contains the scripts used to analyzed >600 complete genomes of *Aedes aegypti* (*Ae. aegypti*) from 40 populations worldwide to uderstand the origin and evolution of its **domestication**. The repository is divided in 4 sections: 1) `processing samples`, 2) `variant calling`, 3) detection of `selection`, and 4) `Supplementary Data`.



### Processing samples

This directory contains three files:
* **`quality_control.md`:** Commands used to performed a quality control in reads (\*.fastq.gz) and alignments (BAM format).
* **`alignments.sh`:** Bash script (slurm jobs cluster system) used to align paired ended reads to the reference genome of *Ae. aegypti*, AaegL5.
* **`recalibration.md`:** Set of command lines used to 1) re-align indels and 2) recalibrate alignments in BAM format.
  

### Variant calling using (GATK)

Single bash script `variant_calling_gatk.sh`. It was used in each population for:
  1) Call variants (indels & SNPs). 
  2) SNPs extraction.
  3) Filter SNPs in two step procedure to get a final set of biallelic SNPs.


### Selection

This directory contains different scripts:\
<ins>**PCadapt:**</ins>
* **`outliers_bestK_pcadapt.R`:** R script to identify the best K used to identify outlier SNPs across populations using PCadapt.
* **`outliers_pcadapt.R`:** R script to identify the total set of outlier of SNPs of a population using PCadapt.

<ins>**RAiSD: FDR of 5%**</ins>
* **`RAiSD_hard_sweeps.md`:** Descriptive steps to do RAiSD analysis at local or global level in our samples.
* **`outliers_RAiSD.R`:** R script to identify the total set of outliers in each RAiSD analysis at local or global level.

<ins>**RAiSD: Top 0.05% and 1% of tail `mu` distribution:**</ins>
* **`selective_sweeps.raisd.LOCAL.get_tail_tops_1_n_0.05.R` :** R script that performs two main cutoffs and creates two BED outfiles to map against the reference genome of *Ae. aegypti* (AaegL5).
* **`gene_distribution_across_Africa_and_Out_of_Africa.LOCAL_RAiSD_analysis.pl` :** PERL script that filters the output files of each population to define cases that in which genesa are either shared and/or specific for a  population. It focused on the cases of Out of Africa specific, but also shows cases for shared between Africa and Out of Africa, and Africa population specific.

<ins>**Detection of Selection using the McDonald and Kreitman test (MKT) and the Direction of Selection (DoS):**</ins>

Summary steps description and code is in the page: `Dection_of_Selection_with_MKT.md`. 


## 
## Supplementary Data: 

This directory contains the files related to **Supplementary Data** associated to the current submitted manuscript.

* Supplementary Data Files 1-7 in proper format and stored in a zipped file: `Supplementary_Data_Files.zip`: 

 1) Supplementary Data 1. SNP statistics for populations through genomic regions (TXT).
 2) Supplementary Data 2. Sequences of new detected nrEVEs (FASTA).
 3) Supplementary Data 3. Phylogenetic trees for populations and individuals (NEWICK).
 4) Supplementary Data 4. Information for 8,148 hard selective sweeps detected with RAiSD in out-of-Africa populations (TXT).
 5) Supplementary Data 5. Information for 1,030 SNP outliers detected with PCAdapt within 2,266 genes (VCF format).
 6) Supplementary Data 6. Matrix with DoS scores for 11,651 orthologous protein-coding genes in AaegL5 and each Ae. aegypti population (TXT).
 7) Supplementary Data 7. Matrix with MKT scores for 11,651 orthologous protein-coding genes in AaegL5 and each Ae. aegypti population (TXT).
 8) Supplementary Data 8. Matrix with DoS scores used to estimate relaxed selection (TXT).
 9) Supplementary Data 9. Genomic coordinates of SNPs in AaegL5 obtained from the literature and VectorBase (TXT).
 10) Supplementary Data 10. Source data of metrics used to plot Figure 4b (TXT).

