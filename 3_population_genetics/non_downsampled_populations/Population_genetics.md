
## SNP counts and other Population Genetics analyses

This file contains three main sections:

* **`SNP density:`** Quantification of the total number of SNPs across the whole genome (chromsomes only) based on different non-overlapping sliding windows using `VCFtools`.
* **`Nucletoide diversity:`** Quantification of the nucleotide diversity (pi) across the whole genome (chromsomes only) based on different non-overlapping sliding windows using `VCFtools`.
* **`Tajima's D:`** Quantification of the Tajima's D across the whole genome (chromsomes only) based on different non-overlapping sliding windows using `VCFtools`.

* **`PCA and Admixture:`** Population structure analysis using the SNPs datasets from all population merged into a single VCF file.


**NOTE:** These series of commands, pipelines and/or scripts show the steps used for our analyses. They were used repeatedly and separately for our (a) different SNPs datasets (e.g. NR-SNPs, R-SNPS, EXONs, whole genome), (b) our different slinding window sizes analyses (500, 250, 100, 50 and 10 kb), (c) for a single population or (d) for population strucutre analysis (all SNPs datasets merged into a single VCF file). Therefore, to avid a high redundancy in all command pipelines and to not make an extense readme file, we show here one single example 'file name' case: the "NR-SNPs" dataset (SNPs dataset without repeated sequences).

----

### SNP density

Single population analysis.

```bash
PATH_to_myVCFs=/full/path/directory/all_populations/VCF_files/;
outfilesDIR=/full/path/directory/all_populations/metrics_outfiles;

# make a list of jobs by listing each population VCF filtered file name
ls $PATH_to_myVCFs/*.fltrd_snps.repeats_removed.biallelic.vcf.gz | sed 's/.vcf.gz//' | parallel --dryrun 'vcftools --gzvcf {}.vcf.gz --SNPdensity 50000 --out $outfilesDIR/{/}.snpDensity.50kb' > aedes_aegypti.all_vcfs_populations.nr_snps.SNPdensity.50kb_jobs.sh;

# run all jobs of each single list in parallel using 40 cores (1 core/pop):
for myJOBS in aedes_aegypti.all_vcfs_populations.nr_snps.SNPdensity.*_jobs.sh; do parallel --jobs 40 < $myJOBS; done
```



### Nucleotide diversity

Single population analysis.

```bash
PATH_to_myVCFs=/full/path/directory/all_populations/VCF_files/;
outfilesDIR=/full/path/directory/all_populations/metrics_outfiles;


# make a list of jobs by listing each population VCF filtered file name
ls $PATH_to_myVCFs/*.fltrd_snps.repeats_removed.biallelic.vcf.gz | sed 's/.vcf.gz//' | parallel --dryrun 'vcftools --gzvcf {}.vcf.gz --window-pi 50000 --window-pi-step 50000 --out $outfilesDIR/{/}.nucdiv.50kb' > aedes_aegypti.all_vcfs_populations.nr_snps.nucdiv.50kb_jobs.sh;

# run all jobs of each single list in parallel using 40 cores (1 core/pop):
for myJOBS in aedes_aegypti.all_vcfs_populations.nr_snps.nucdiv.*_jobs.sh; do parallel --jobs 40 < $myJOBS; done
```



### Tajima's D

Single population analysis.

```bash
PATH_to_myVCFs=/full/path/directory/all_populations/VCF_files/;
outfilesDIR=/full/path/directory/all_populations/metrics_outfiles;


# make a list of jobs by listing each population VCF filtered file name
ls $PATH_to_myVCFs/*.fltrd_snps.repeats_removed.biallelic.vcf.gz | sed 's/.vcf.gz//' | parallel --dryrun 'vcftools --gzvcf {}.vcf.gz --TajimaD 50000 --out $outfilesDIR/{/}.tajsD.50kb' > aedes_aegypti.all_vcfs_populations.nr_snps.tajsD.50kb_jobs.sh;

# run all jobs of each single list in parallel using 40 cores (1 core/pop):
for myJOBS in aedes_aegypti.all_vcfs_populations.nr_snps.tajsD.*_jobs.sh; do parallel --jobs 40 < $myJOBS; done
```

----

### Principal Component Analysis

Global populations analysis (single VCF file where all populations were merged). Coverage of >90% individuals per population.

```bash
# [1] After merging all VCF files with BCFtools (merge function) and 
# [2] retain all SNPs present in at least 90% of all populations using PLINK, what is next is [3] to use the function "pca" of PLINK to create the main outfiles that are required to make the PCA plot: 'eigenvec' and 'eigenval'.

# 1. Merge VCFs (a specific set of genomic data, such as "SNPs without repeated sequences" (a.k.a NR-SNPs, as named in the paper):
ls $PATH_to_myVCFs/*.fltrd_snps.repeats_removed.biallelic.vcf.gz > aedes_aegypti.all_pops.nr_snps.biallelic.vcfs_list.txt;

# this will take a long time (but it is a way faster (2x-to-3x) than other programs!!)
bcftools merge -l aedes_aegypti.all_pops.nr_snps.biallelic.vcfs_list.txt -Oz -o aedes_aegypti.all_pops_merged.nr_snps.biallelic.vcf.gz;

bcftools index -t aedes_aegypti.all_pops_merged.nr_snps.biallelic.vcf.gz;


# 2. Identify all linked alleles (Linkage Disequilibrium, LD), options:
plink \
  --vcf aedes_aegypti.all_pops_merged.nr_snps.biallelic.vcf.gz \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --indep-pairwise 50 10 0.1  --make-bed --max-alleles 2  --geno 0.1 \
  --threads 30 \
  --out $outPCAdir/aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01


# 2. Thin the data by filtering out linked alleles (LD), options:
plink \
  --vcf aedes_aegypti.all_pops_merged.nr_snps.biallelic.vcf.gz \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --extract $outPCAdir/aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.prune.in \
  --threads 30 \
  --make-bed  --pca --out $outPCAdir/aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned
wait;


# Final notes:
# It will generate two main files:
#   eigenvec: aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.eigenvec
#   eigenval: aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.eigenval

# Use both outfiles to make the PCA plot.

# Repeat this procedure for each genomic region to be analyzed: repeated sequences, non-repeated sequences, exons, whole genome, etc.
```



### Admixture

Global populations analysis (single VCF file where all populations were merged). Coverage of >90% individuals per population.

```bash
# We used the program 'Admixture' and we transformed the VCF (same used for the PCA) to a PLINK format. Then, we ran the program for each of our SNPs datasets separately (e.g. NR-SNPs, R-SNPS, EXONs, whole genome). Next, we evaluated the cross-validation error for each of these runs. To not make repetitive all command pipelines one single example 'file name' case is used: the "NR-SNPs" dataset.


# First: "reformat" the infile *.bim file (from plink) to be read it by admixture
# ---------------------------------------------------------------------------------------------
# Why? ADMIXTURE does not accept chromosome names that are not in a format of human chromosomes: we just need to exchange the first column by 0

# dorectory where all admixture's outfiles will be stored. 
admixtureDIR=`pwd`;

# taking infile from other directory
awk '{$1=0;print $0}' $outPCAdir/aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.bim  > $admixtureDIR/aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.bim.tmp;
# rename
cd $admixtureDIR;
mv aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.bim.tmp   aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.bim;
# only read
chmod 0444 aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.bim;

cp $outPCAdir/aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.bed .


# Second: make admixture jobs
# ---------------------------------------------------------------------------------------------
# syntaxis command with range of 40 clusters (k=1..40), job-specific SEED (time linux command), 100 bootstraps, and 40 cores.
#   admixture --cv -B1000 -j20 -s time snps.vcf2plink.pruned.bed 2>> snps.vcf2plink.pruned.admixture.boots.kxxx.log > snps.vcf2plink.pruned.admixture.boots.kxxx.log
# where: xxx = k number

# jobs list (each with 40 cores)
for K in $(seq 40); do echo "admixture --cv -B1000 -j40 -s time aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.bed $K > aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.pruned.k$K.log"; done > snps.biallelic.geno90.LD-01.prune_cleaned.ADMIXTURE.jobs.sh

# Run a single job in the jobs list:
parallel --joblog snps.biallelic.geno90.LD-01.prune_cleaned.ADMIXTURE.jobs.runtime.log  --jobs 1 < snps.biallelic.geno90.LD-01.prune_cleaned.ADMIXTURE.jobs.sh 2>> snps.biallelic.geno90.LD-01.prune_cleaned.ADMIXTURE.jobs.stderr.log >> snps.biallelic.geno90.LD-01.prune_cleaned.ADMIXTURE.jobs.stderr.log &


# Third: get cross-validation (CV)
# ---------------------------------------------------------------------------------------------
# get CV error per ran K:
grep -h CV $admixtureDIR/*.log > aedes_aegypti.all_pops_merged.nr_snps.biallelic.geno90.LD-01.prune_cleaned.ADMIXTURE.CV_scores.txt

# read $admixtureDIR/*.Q files to make the admixture stacked-bar plots
```

