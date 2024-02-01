
## Calculate the genetic diversity and Tajima's D

**1) Calculate genotype likelihoods and the site allele frequency (SAF)** 

```bash
myWDIR=/scr/core/nlozada/aedes_aegypti/scripts/;
outDIR=/scr/core/nlozada/aedes_aegypti/output;
popLIST=my_aaeg_populations.list.txt;   # list with a full path for a TXT file containing the WGS BAM alignments of each inidividual in a population 
refGenome=/scr/core/nlozada/aedes_aegypti/reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta;   # reference genome (in fasta format) 


# calculate genotype likelihoods using the GATK model

# main output extension will be *.saf
for i in $popLIST; do
   POP=$(echo $i | cut -d'.' -f 1);  
   angsd -b $i -anc $refGenome -out $outDIR/${POP}.angsd -ref $refGenome -minMapQ 10 -minQ 10 -minInd 1 -doSaf 1 -GL 2 -nThreads 8;
   wait;
   sleep 2;
done

wait;
```

**2) Calculate the site frequency spectrum (SFS)**

```bash
popSAFindex=$outDIR/*.angsd.saf.idx;

# main output extension will be *.sfs
for i in $popSAFindex; do
   POP=$(echo $i | cut -d'.' -f 1);
   realSFS  $i  -maxiter 100  -cores 60  >  ${POP}.saf2sfs.sfs 2>&1 | tee ${POP}.saf2sfs.stderr.log;
   wait;
   sleep 2;
done

wait;
```
  
**3) Calculate the thetas (population scaled mutation rate) for each site**

Based on these mutation rates for each population, genetic diversity and divergence metrics can be calculated using subprograms `saf2theta` and `do_stat` based on [Korneliussen et al., 2013](https://doi.org/10.1186/1471-2105-14-289).


```bash
popSAFlist=scr/core/nlozada/aedes_aegypti/outputs/*.angsd.saf.idx;
popSFSlist=scr/core/nlozada/aedes_aegypti/outputs/*.angsd.saf2sfs.sfs;

# get population name from *.saf files
for i in $popSAFlist; do
    POPa=$(echo $popSAFlist | cut -d'.' -f 1);

    # get population name from *.sfs files
    for j in $popSFSlist; do
        POPb=$(echo $popSFSlist | cut -d'.' -f 1);

        if [[ "$POPa" == "$POPb" ]] {
           realSFS saf2theta $i  -sfs $j  -outname $outDIR/${POPb}.thetas.stdout.txt 2>>  $outDIR/${POPb}.thetas.stderr.log;
           wait;
           sleep 2;
        }
    done;
    # global statistics:
    thetaStat do_stat  $outDIR/${POPa}.thetas.stdout.txt > $outDIR/${POPa}.thetas.global_stats.stdout.txt 2>>  $outDIR/${POPa}.thetas.global_stats.stderr.log;
    # theta site specific:
    thetaStat print  $outDIR/${POPa}.thetas.stdout.txt > $outDIR/${POPa}.thetas.persite_stats.stdout.txt 2>>  $outDIR/${POPa}.thetas.persite_stats.stderr.log;
done

wait;
```

**4) Calculate summary statistics per pupulation**

Global theta summary basic statitics per population were summarize per `nuncleotide diversity (π)` and `Tajima's D` using a custom R script that simply used two main functions, `group_by` and `summarize` from the [R `dplyr` package version 1.1.4](https://dplyr.tidyverse.org).

```R
library(dplyr)

# -------------------------------------------------------------------------------------------------
infile_thetas_table <- read.csv("thetas.stats.global.txt", sep="\t", headers=TRUE);


alpha=0.05;

# nucleotide diversity (π)
# -------------------------------------------------------------------------------------------------
nucdiv_stats.population_chrms_thetas_scores <- as.data.frame(infile_thetas_table %>% 
                                                                    group_by(POPULATION) %>% 
                                                                               summarise(min  = min(NUC_DIVERSITY),
                                                                                         max  = max(NUC_DIVERSITY),
                                                                                         mean = mean(NUC_DIVERSITY),
                                                                                         sd   = sd(NUC_DIVERSITY),
                                                                                         se   = sd(NUC_DIVERSITY)/sqrt(length(NUC_DIVERSITY)),
                                                                                         q1   = quantile(NUC_DIVERSITY, 0.25),
                                                                                         q3   = quantile(NUC_DIVERSITY, 0.75),
                                                                                         t.score = qt(p=alpha/2, df=df,lower.tail=F),
                                                                                         margin.error = t.score * se,
                                                                                         ci.lower = mean - margin.error,
                                                                                         ci.upper = mean + margin.error
                                                                                ) 
                                             )

# my stats :
nucdiv_stats.population_chrms_thetas_scores



# Tajima's D
# -------------------------------------------------------------------------------------------------
tajsD_stats.population_chrms_thetas_scores <- as.data.frame(infile_thetas_table %>% 
                                                                    group_by(POPULATION) %>% 
                                                                               summarise(min  = min(TAJIMASD),
                                                                                         max  = max(TAJIMASD),
                                                                                         mean = mean(TAJIMASD),
                                                                                         sd   = sd(TAJIMASD),
                                                                                         se   = sd(TAJIMASD)/sqrt(length(TAJIMASD)),
                                                                                         q1   = quantile(TAJIMASD, 0.25),
                                                                                         q3   = quantile(TAJIMASD, 0.75),
                                                                                         t.score = qt(p=alpha/2, df=df,lower.tail=F),
                                                                                         margin.error = t.score * se,
                                                                                         ci.lower = mean - margin.error,
                                                                                         ci.upper = mean + margin.error
                                                                                ) 
                                             )

# my stats :
tajsD_stats.population_chrms_thetas_scores

```
