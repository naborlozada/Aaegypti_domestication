# Aaegypti_domestication
---
Title:
**Genomic signatures of globally invasive *Aedes aegypti* populations***

**Authors:** Alejandro Nabor Lozada-Chávez, Irma Lozada-Chávez, Niccolò Alfano, Umberto Palatini, Davide Sogliani, Samia Elfekih, Teshome Degefa, Maria V. Sharakhova, Athanase Badolo, Sriwichai Patchara, Mauricio Casas-Martinez, Bianca C. Carlos, Rebeca Carballar-Lejarazú, Luis Lambrechts, Jayme A. Souza-Neto & Mariangela Bonizzoni

**CURRENTLY, THIS MANUSCRIPT IS UNDER REVIEW.**

--- 

<ins>These and previous versions of the scripts are subject to changes and tests. If you use them (one or more) in a complete or partial form for your research, as well as any other data contained here, please cite this repository (**`Github` and `Zenodo`**) and the preprint **bioRxiv** paper.</ins>

  **Repository/scripts citation:**\
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10721413.svg)](https://doi.org/10.5281/zenodo.7863455)
   
   * Nabor Lozada. (2024). Aedes aegypti domestication (v0.4). Zenodo. https://doi.org/10.5281/zenodo.10721413
   * A. N. Lozada-Chávez et. al. 2023. [Molecular signature of domestication in the arboviral vector *Aedes aegypti*](https://doi.org/10.1101/2023.03.13.532092). bioRxiv. DOI: https://doi.org/10.1101/2023.03.13.532092 [submitted]. 


##
## NOTICE (updated 3-October-2024):
1) **selection** In this directory, a single MD file describing the approach to detect Selection on protein coding genes based on both, the McDonald and Kreitman test (MKT) and the Direction of Selection (DoS). First, single short text steps are decribed from the detection of orthologs, codon alignments, and both analyses of MKT and DoS. Next, a summary of each step in bash scripting languages for orthologs detection and codon alignments, and both analyses, MKT and DoS, in R code context MD format.
2) **Supplementary Data** associated to the current submmitted version of the manuscript was updated and is located in the <ins>**Supplementary Data** directory</ins> of this repository. A summary of the content is at the bottom of this page. 
##
## NOTICE (updated 13-May-2024):
1) **selection** In this directory, six scripts and a summary step protocol were added, and they are related to the <ins>**detection of outliers based on the selection of a top 0.05% and 1% of a tail distribution** of `RAiSD` results, and to the **detection of selection in coding protein genes using PAML**. See **selection** directory for details.</ins>.
2) The Supplementary Data ZIP file was updated. See changes in the **Supplementary Data** directory.
3) SRA name of each samples were added to the Supplementary Table 1 (BioProject: PRJNA943178).
##
## NOTICE (updated 15-Jan-2024):
1) **Downsampling** description files and scripts were included in this version in the <ins>**downsampled_analyses** directory</ins>.
2) **Population Branch Statistics  (PBS)** description file and scripts were included in this version in the <ins>**downsampled_analyses** directory</ins>.
3) **Detection of hard selective sweeps** description files and scripts associated to the new integrated approach **RAiSD** were included in the <ins>**Selection** directory</ins>. 
##
## NOTICE (updated 22-May-2023):
1) **Supplementary Data** associated to the current submmitted version of the manuscript is in the <ins>**Supplementary Data** directory</ins> of this repository. A summary of the content is below. 
2) Due technical issues, we are currently reuploading the Whole Genome Sequences (WGS) in the SRA database under the BioProject: PRJNA943178.

## 
## Description

Scripts used to analyzed >600 complete genomes of *Aedes aegypti* (*Ae. aegypti*) from 39 populations worldwide to uderstand the origin and evolution of its **domestication**. Sections of this repository: `processing samples`, `variant calling`, detection of `selection`.



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

<ins>**pN/ps ratio:**</ins>\
(1) Based on only segregating sites of Ae. aegypti's populations:
* **`pnps_ratio.R`:** R script to calculate the pn/pn ratio of SNPs of a population.

(2) Based on the divergence between two sister species: Ae. aegypti and Ae. albopictus using YN model from PAML. 
* **`detection_of_selection_PAML.md`** Summary protocol steps and the script used to calculate the pN/pS ratio for each protein coding gene in each single population.
* **`merge_orthologues_genes.pl`** PERL script ot join two orthologs genes in a single fasta file.
* **`GENEID_yn00.ctl`** Control infile from PAML YN model (YN_model.ctl). 
* **`get_pnps_ratio_per_poplation.pl`** PERL script that extract the pN/pS statistics from PAML analysis and sort them for each single gene present in a population, and define the type of selection based on the threshold.


### Other scripts

This directory contains scripts used for additional analyses and/or parse several results. 


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

