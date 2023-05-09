# Aaegypti_domestication

Title:
**Molecular signature of domestication in the arboviral vector *Aedes aegypti***

**Authors:** Alejandro Nabor Lozada-Chávez, Irma Lozada-Chávez, Niccolò Alfano, Umberto Palatini, Davide Sogliani, Samia Elfekih, Teshome Degefa, Maria V. Sharakhova, Athanase Badolo, Sriwichai Patchara, Mauricio Casas-Martinez, Bianca C. Carlos, Rebeca Carballar-Lejarazú, Luis Lambrechts, Jayme A. Souza-Neto & Mariangela Bonizzoni


<ins>These and previous versions of the scripts are subject to changes and tests. If you use them (one or more) in a complete or partial form for your research, as well as any other data contained here, please cite this repository (**`Zenodo`**) and the preprint **bioRxiv** paper.</ins>

  **Repository/scripts citation:**\
   [![DOI](https://zenodo.org/badge/630563603.svg)](https://zenodo.org/badge/latestdoi/630563603)
   * Nabor Lozada. (2023). naborlozada/Aaegypti_domestication: Aaegypti_domestication (v0.1b). Zenodo. https://doi.org/10.5281/zenodo.7863456
   * A. N. Lozada-Chávez et. al. 2023. [Molecular signature of domestication in the arboviral vector *Aedes aegypti*](https://doi.org/10.1101/2023.03.13.532092). bioRxiv. DOI: https://doi.org/10.1101/2023.03.13.532092 [submitted]. 

##
## NOTICE:
Due technical issues, we are currently reuploading the Whole Genome Sequences (WGS) in the SRA database under the BioProject: PRJNA943178.  


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




## 
## Supplementary Information: 

This directory contains the files related to **Supplementary Information** and **Supplementary Data Files** associated to the current submitted manuscript.

* `Supplementary_Information.pdf`: Information merged into a single PDF: Supplementary Information, Supplementary Data 4, and Supplementary Data 6.
* `Supplementary_Data_Files.zip`: Supplementary Data Files 1-6 in proper format.


