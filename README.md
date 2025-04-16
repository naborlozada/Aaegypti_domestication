# Aaegypti_domestication
---
Title:
**Adaptive genomic signatures of globally invasive popultions of the yellow fever mosquito *Aedes aegypti****


\* These authors contributed equally in the research article: Alejandro N. Lozada-Chávez, Irma Lozada-Chávez.


**Citation:**\
Lozada-Chávez AN*, Lozada-Chávez I*, Alfano N, Palatini U, Sogliani D, Elfekih S, Degefa T, Sharakhova MV, Badolo A, Sriwichai P, Casas-Martínez M, Carlos BC, Carballar-Lejarazú R, Lambrechts L, Souza-Neto JA, Bonizzoni M. 2025. Adaptive genomic signatures of globally invasive populations of the yellow fever mosquito *Aedes aegypti*. ***Nat. Ecol. Evol.*** 9(4):652-671. doi: [10.1038/s41559-025-02643-5.](https://www.nature.com/articles/s41559-025-02643-5)

PMID: 40155778; PMCID: PMC11976285. Early publication: 2025 Mar 28.



--- 
<br>
<br>

* This final version of the repository has the main scripts used for our study (See Description). 
* All Supplementary and Additional Information, such as Tables (single excel file), Extended Data Figures and more can be found in the *Supplementary Information* section of our paper in [**Nature Ecology & Evolution**](https://www.nature.com/articles/s41559-025-02643-5#Sec29).
* Supplementary Data (1-12) is accessible from the Zenodo repository. The final and complete version of these datasets is in the **version 3**, and is accessible from this DOI:\
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14948092.svg)](https://doi.org/10.5281/zenodo.14948092)
* **NOTICE:** *On April 2nd, 2025*, we performed a slight reorganization of this repository by separating in two different directories the `Populations Genetics` and the `Detection of Selection` sections. Also, we arranged  properly few file names used in some commands to make consistent the steps followed one-after-another.   

<br>

**IMPORTANT NOTE:** If you use one or more of these scripts and/or data in a complete or partial form for your research, as well as any other information contained here, **please cite us with at least one of these references based on:**

   * **Research paper:** Lozada-Chávez, A.N., Lozada-Chávez, I. *et al.* Adaptive genomic signatures of globally invasive populations of the yellow fever mosquito Aedes aegypti. Nat Ecol Evol (2025). DOI: https://doi.org/10.1038/s41559-025-02643-5
   * **Scripts and data in Github**:  Lozada-Chávez, A. N. 2025. Aaegypti_domestication. Github repository. https://github.com/naborlozada/Aaegypti_domestication.

   * **Datasets:** Lozada-Chávez, A. N. et al. Adaptive genomic signatures of globally invasive populations of the yellow fever mosquito *Aedes aegypti*. Zenodo https://doi.org/10.5281/zenodo.14948092 (2024).
   

<br>
For historical reasons, here we cite our preprint in **bioRxiv** (not peer reviewed!):

* A. N. Lozada-Chávez et. al. 2023. Molecular signature of domestication in the arboviral vector *Aedes aegypti*. bioRxiv. DOI: https://doi.org/10.1101/2023.03.13.532092. 


## 
## Description

This repository contains the scripts used to analyzed 686 complete genomes of *Aedes aegypti* (*Ae. aegypti*) from 40 populations worldwide to uderstand the origin and evolution of its **domestication**. The repository is divided in 4 sections: 1) `processing samples`, 2) `variant calling`, 3) `population genetics`, and 4) `detection of selection`.



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


### Population Genetics analyses

  1) Analyses using `non-downsampled populations`: SNP density, Nucletide diversity, Tajima's D, PCA, and Admixture.  
  2) Analyses based on the `downsampled populations`: Nucletide diversity, Tajima's D, and Population Branch Statics (PBS) using ANGSD.

NOTE: In the section `downsampled populations` we provided the full original output data from ANGSD for nucleotide diversity and Tajima's D by population, as well as a simple R script to produce basic statistics, so the results in Table 1 and Supplementary Table 5 (tables 2 and 3) from our paper can be reproduced. The original data from ANGSD is also shared in the Supplementary Table 6, however, we noticed that data of these two metrics (nucleotide diversity and Tajima's D) are missing for the population Entebbe. Now are included in the TXT file attached in this section).

### Detection of Selection

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

This **Supplementary Data** can be found in **ZENODO**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14948092.svg)](https://doi.org/10.5281/zenodo.14948092)


Supplementary_Data.zip: Supplementary Data Files 1-12 in compressed ZIP or TAR formats. The Supplementary Data 12 (SD-12) is missing in the Supplementary Data in Nature Ecology and Evolution since it was a large big file. In Zenodo SD-12 can be found. A complete list of SDs is below:

1) Supplementary Data 1. SNP statistics for populations through genomic regions (TXT). 
2) Supplementary Data 2. Sequences of new detected nrEVEs (FASTA). 
3) Supplementary Data 3. Phylogenetic trees for populations and individuals (NEWICK). 
4) Supplementary Data 4. Information for 8,120 hard selective sweeps detected with RAiSD in out-of-Africa populations (TXT). 
5) Supplementary Data 5. Information for 1,030 SNP outliers detected with PCAdapt within 2,266 genes (VCF format). 
6) Supplementary Data 6. Matrix with DoS scores for 11,651 orthologous protein-coding genes in AaegL5 and each Ae. aegypti population (TXT). 
7) Supplementary Data 7. Matrix with MKT scores for 11,651 orthologous protein-coding genes in AaegL5 and each Ae. aegypti population (TXT). 
8) Supplementary Data 8. Matrix with DoS scores used to estimate relaxed selection (TXT). 
9) Supplementary Data 9. Matrix with SNPs and genomic coordinates within adaptive protein-coding genes and ncRNAs that are shared or private for out-of-Africa populations against African populations (TXT). 
10) Supplementary Data 10. Matrix with 483 nonsynonymous SNPs and their allele frequencies for our 40 populations Florida and Colombia (TXT).
11) Supplementary Data 11. Genomic coordinates of SNPs in AaegL5 obtained from the literature and VectorBase (TXT). 
12) Supplementary Data 12. Source data of metrics used to plot Figure 4b (TXT).


NOTE: See notes added on the `description` in the Zenodo's repository.

