# Aaegypti_domestication
---
Title:
**Genomic signatures of globally invasive *Aedes aegypti* populations***

**Authors:** Alejandro Nabor Lozada-Chávez, Irma Lozada-Chávez, Niccolò Alfano, Umberto Palatini, Davide Sogliani, Samia Elfekih, Teshome Degefa, Maria V. Sharakhova, Athanase Badolo, Sriwichai Patchara, Mauricio Casas-Martinez, Bianca C. Carlos, Rebeca Carballar-Lejarazú, Luis Lambrechts, Jayme A. Souza-Neto & Mariangela Bonizzoni

**CURRENTLY, THIS MANUSCRIPT IS UNDER REVIEW.**

--- 

<ins>These and previous versions of the scripts are subject to changes and tests. If you use them (one or more) in a complete or partial form for your research, as well as any other data contained here, please cite this repository (**`Github` and `Zenodo`**) and the preprint **bioRxiv** paper.</ins>

  **Repository/scripts citation:**\
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10721413.svg)](https://doi.org/10.5281/zenodo.10721413)
   
   * Nabor Lozada. (2024). Aedes aegypti domestication (v0.4). Zenodo. https://doi.org/10.5281/zenodo.10721413
   * A. N. Lozada-Chávez et. al. 2023. [Molecular signature of domestication in the arboviral vector *Aedes aegypti*](https://doi.org/10.1101/2023.03.13.532092). bioRxiv. DOI: https://doi.org/10.1101/2023.03.13.532092 [submitted]. 

##
## NOTICE (updated 13-May-2024):
1) **selection** In this directory, two scripts were added, and they are related to the <ins>**detection of outliers based on the selection of a top 0.05% and 1% of a tail distribution** of `RAiSD` results</ins>.
     * `R` script: performs the cutoffs and mapping.
     * `PERL` script. Analyzed the gene distribution, specific and/or shared, among Out of Africa and Africa populations.
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

This directory contains two files:
* **`outliers_bestK_pcadapt.R`:** R script to identify the best K used to identify outlier SNPs across populations using PCadapt.
* **`outliers_pcadapt.R`:** R script to identify the total set of outlier of SNPs of a population using PCadapt.
* **`pnps_ratio.R`:** R script to calculate the pn/ps ratio of SNPs of a population.
* **`RAiSD_hard_sweeps.md`:** Descriptive steps to do RAiSD analysis at local or global level in our samples.
* **`outliers_RAiSD.R`:** R script to identify the total set of outliers in each RAiSD analysis at local or global level.

### Other scripts

This directory contains scripts used for additional analyses and/or parse several results. 


## 
## Supplementary Data: 

This directory contains the files related to **Supplementary Data** associated to the current submitted manuscript.

* Supplementary Data Files 1-7 in proper format and stored in a zipped file: `Supplementary_Data_Files.zip`: 

 1) Supplementary Data 1: Population information and basic SNP statistics (TXT file).
 2) Supplementary Data 2: new nrEVEs in `fasta` sequences (TXT file).
 3) Supplementary Data 3: two phylogenetic trees, based on individuals (RAxML) and populations (Treemix), are provided in `newick` format (TXT file).
 4) Supplementary Data 4: 10,030 SNP outliers associated to 2,266 genes are provided in `VCF` format (\*.VCF.GZ file).
 5) Supplementary Data 5: Matrix with the pN/pS gene ratio across the whole genome of *Ae. aegypti* using the improved method of Li's Kimura 2 Parameters (K2P) approach (TXT file).
 6) Supplementary Data 6: SNP dataset obtained from the literature to recalibrate genomes (TXT file).
 7) Supplementary Data 7: Table with the pN/pS gene ratio across the whole genome of *Ae. aegypti* using the PAML protocol approach and the *Ae. albopictus* species as a sister species (TXT file).

