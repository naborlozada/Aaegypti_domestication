
## Calculate the genetic diversity and Tajima's D in the genome of *Aedes aegypti* for all our populations

For each population, the site allele frequency (SAF), site frequency spectrum (SFS), and the thetas (population scaled mutation rate) for each site were calculated using the command lines shown below. All programs executions (runs below) were implemented in parallel mode in different servers and/or clusters with high computational performance to accelerate the calculations.

**1) Calculate genotype likelihoods and the site allele frequency (SAF)** 

```bash
myWDIR=/scr/core/nlozada/aedes_aegypti/scripts/;
outDIR=/scr/core/nlozada/aedes_aegypti/output;
popLIST=my_aaeg_populations.list.txt;              # list with a full path for a TXT file containing the WGS BAM alignments file names of each inidividual in a population 
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

# read SAF file names
for i in $popSAFlist; do
    #get population name from each *.saf file
    POPa=$(echo $popSAFlist | cut -d'.' -f 1);

    # read SFS file names
    for j in $popSFSlist; do
        # get population name from each *.sfs file
        POPb=$(echo $popSFSlist | cut -d'.' -f 1);
        
        if [[ "$POPa" == "$POPb" ]]; then
           realSFS saf2theta $i  -sfs $j  -outname $outDIR/${POPb}.thetas.stdout.txt 2>>  $outDIR/${POPb}.thetas.stderr.log;
           wait;
           sleep 2;
        fi
    done;
    # global statistics:
    thetaStat do_stat  $outDIR/${POPa}.thetas.stdout.txt > $outDIR/${POPa}.thetas.global_stats.stdout.txt 2>>  $outDIR/${POPa}.thetas.global_stats.stderr.log;
    # theta site specific:
    thetaStat print  $outDIR/${POPa}.thetas.stdout.txt > $outDIR/${POPa}.thetas.persite_stats.stdout.txt 2>>  $outDIR/${POPa}.thetas.persite_stats.stderr.log;
done

wait;
```

**4) Calculate summary statistics per population and by groups (Africa vs Out of Africa)**

Global theta summary basic statistics per population were summarize for both `nuncleotide diversity (π)` and `Tajima's D` separately using a custom R script that simply uses two main functions, `group_by` and `summarize` from the [R `dplyr` package version 1.1.4](https://dplyr.tidyverse.org), as well as additional parsing functions (select, filter, mutate) in pipeline commands to get (i.e. *reproduce*) the information in the Table 1 and ST-5. The script below uses a TXT file named "aedes_aegypti.downsampled_populations.angsd_pi_tajD_scores_byChrm.parsed.txt", which is attached in this directory and included in the ST-6 from our `Supplementary Tables` excel file. Basic information can be found in the TXT file.

**IMPORTANT NOTE: In the original ST-6 (online Supplmentary Tables), we noticed that the scores of both metrics for the "Entebbe" (ENT) population were missing, so we provided them in this TXT file.**

```R

.libPaths("/full/path/to/R_library/myPackages/");

library(dplyr);


# working directory
workdir <- "/full/path/of/files/";
setwd(workdir);


# cleaning...
rm(list=ls());


# infile
# -------------------------------------------------------------------------------------------------
infile_divmetrics <- read.table(file="aedes_aegypti.downsampled_populations.angsd_pi_tajD_scores_byChrm.parsed.txt", sep="\t", header=TRUE, comment.char = "#");



# ------------------------------------------------------------------------------------------------- #
#                       *** Reproduce stats in ST-5: tables 2 and 3 ***                             #
# ------------------------------------------------------------------------------------------------- #

alpha=0.05;

# nucleotide diversity (π)
# -------------------------------------------------------------------------------------------------
nucdiv_stats.Populations <- infile_divmetrics %>% 
                                 group_by(Population) %>% 
                                       summarise(df = length(Pi_ratio)-1,
                                                   min  = min(Pi_ratio),
                                                   max  = max(Pi_ratio),
                                                   mean = mean(Pi_ratio),
                                                   sd   = sd(Pi_ratio),
                                                   se   = sd(Pi_ratio)/sqrt(length(Pi_ratio)),
                                                   t.score = qt(p=alpha/2, df=df,lower.tail=F),
                                                   margin.error = t.score * se,
                                                   ci.lower = mean - margin.error,
                                                   ci.upper = mean + margin.error
                                          ) %>% as.data.frame()




# Tajima's D
# -------------------------------------------------------------------------------------------------
tajsD_stats.Populations <- infile_divmetrics %>% 
                                          group_by(Population) %>% 
                                                summarise(df = length(Tajima)-1,
                                                            min  = min(Tajima),
                                                            max  = max(Tajima),
                                                            mean = mean(Tajima),
                                                            sd   = sd(Tajima),
                                                            se   = sd(Tajima)/sqrt(length(Tajima)),
                                                            t.score = qt(p=alpha/2, df=df,lower.tail=F),
                                                            margin.error = t.score * se,
                                                            ci.lower = mean - margin.error,
                                                            ci.upper = mean + margin.error
                                                   ) %>% as.data.frame()


# defined datasets
nucdiv_stats.Populations$dataset <- "NucDiv";
tajsD_stats.Populations$dataset  <- "TajD";


# sort populations
populations_sorted <- c("Bantata","Kedougou","Mindin","PK10","Dori","Ouagadougou","Ouahigouya","Larabanga","BoabengFiema","Kintampo","Kumasi","Awka","Benoue","Franceville","Libreville","LopeVillage","Bundibugyo","Entebbe","Karenga","Kichwamba","Virhembe","Kakamega","Mbarakani_village","Ganda","Arabuko","KayaBomu","Kwale","ShimbaHills","Rabai_selv","Thies","Ngoye","Rabai_dom","Tapachula","Bebedouro","Santarem","Jeddah","Samut_Sakhon","Bangkok","Tafuna_Village","Zac_Panda");

# sort column rows (population names) based on the "populations_sorted" list: accordingly to ST-5
nucdiv_stats.Populations_strd  <- nucdiv_stats.Populations[match(populations_sorted, nucdiv_stats.Populations$Population),];
tajsD_stats.Populations_strd   <- tajsD_stats.Populations[match(populations_sorted, tajsD_stats.Populations$Population),];



# /// create groups for further stats calculations ///

# Table 2 (ST-5): nucleotide diversity
nucdiv_stats.Population_groups <- nucdiv_stats.Populations_strd %>% 
                                       dplyr::mutate(POPGROUP = case_when(Population == "Arabuko" ~ "Africa",
                                                                          Population == "Awka" ~ "Africa",
                                                                          Population == "Bantata" ~ "Africa",
                                                                          Population == "Dori" ~ "Africa",
                                                                          Population == "Kakamega" ~ "Africa",
                                                                          Population == "Kichwamba" ~ "Africa",
                                                                          Population == "Larabanga" ~ "Africa",
                                                                          Population == "Mindin" ~ "Africa",
                                                                          Population == "PK10" ~ "Africa",
                                                                          Population == "Thies" ~ "Africa",
                                                                          Population == "Benoue" ~ "Africa",
                                                                          Population == "Franceville" ~ "Africa",
                                                                          Population == "Karenga" ~ "Africa",
                                                                          Population == "Kintampo" ~ "Africa",
                                                                          Population == "Libreville" ~ "Africa",
                                                                          Population == "Ngoye" ~ "Africa",
                                                                          Population == "Rabai_dom" ~ "Africa",
                                                                          Population == "ShimbaHills" ~ "Africa",
                                                                          Population == "Virhembe" ~ "Africa",
                                                                          Population == "BoabengFiema" ~ "Africa",
                                                                          Population == "Ganda" ~ "Africa",
                                                                          Population == "KayaBomu" ~ "Africa",
                                                                          Population == "Kumasi" ~ "Africa",
                                                                          Population == "LopeVillage" ~ "Africa",
                                                                          Population == "Ouagadougou" ~ "Africa",
                                                                          Population == "Rabai_selv" ~ "Africa",
                                                                          Population == "Bundibugyo" ~ "Africa",
                                                                          Population == "Kedougou" ~ "Africa",
                                                                          Population == "Kwale" ~ "Africa",
                                                                          Population == "Entebbe" ~ "Africa",
                                                                          Population == "Mbarakani_village" ~ "Africa",
                                                                          Population == "Ouahigouya" ~ "Africa",
                                                                          Population == "Santarem" ~ "Out_of_Africa",
                                                                          Population == "Bangkok" ~ "Out_of_Africa",
                                                                          Population == "Bebedouro" ~ "Out_of_Africa",
                                                                          Population == "Tafuna_Village" ~ "Out_of_Africa",
                                                                          Population == "Zac_Panda" ~ "Out_of_Africa",
                                                                          Population == "Jeddah" ~ "Out_of_Africa",
                                                                          Population == "Samut_Sakhon" ~ "Out_of_Africa",
                                                                          Population == "Tapachula" ~ "Out_of_Africa",
                                                                          TRUE  ~ "warning")
                                          ) %>% as.data.frame();



# Table 3 (ST-5): Tajima's D
tajsD_stats.Population_groups <- tajsD_stats.Populations_strd %>% 
                                       dplyr::mutate(POPGROUP = case_when(Population == "Arabuko" ~ "Africa",
                                                                          Population == "Awka" ~ "Africa",
                                                                          Population == "Bantata" ~ "Africa",
                                                                          Population == "Dori" ~ "Africa",
                                                                          Population == "Kakamega" ~ "Africa",
                                                                          Population == "Kichwamba" ~ "Africa",
                                                                          Population == "Larabanga" ~ "Africa",
                                                                          Population == "Mindin" ~ "Africa",
                                                                          Population == "PK10" ~ "Africa",
                                                                          Population == "Thies" ~ "Africa",
                                                                          Population == "Benoue" ~ "Africa",
                                                                          Population == "Franceville" ~ "Africa",
                                                                          Population == "Karenga" ~ "Africa",
                                                                          Population == "Kintampo" ~ "Africa",
                                                                          Population == "Libreville" ~ "Africa",
                                                                          Population == "Ngoye" ~ "Africa",
                                                                          Population == "Rabai_dom" ~ "Africa",
                                                                          Population == "ShimbaHills" ~ "Africa",
                                                                          Population == "Virhembe" ~ "Africa",
                                                                          Population == "BoabengFiema" ~ "Africa",
                                                                          Population == "Ganda" ~ "Africa",
                                                                          Population == "KayaBomu" ~ "Africa",
                                                                          Population == "Kumasi" ~ "Africa",
                                                                          Population == "LopeVillage" ~ "Africa",
                                                                          Population == "Ouagadougou" ~ "Africa",
                                                                          Population == "Rabai_selv" ~ "Africa",
                                                                          Population == "Bundibugyo" ~ "Africa",
                                                                          Population == "Kedougou" ~ "Africa",
                                                                          Population == "Kwale" ~ "Africa",
                                                                          Population == "Entebbe" ~ "Africa",
                                                                          Population == "Mbarakani_village" ~ "Africa",
                                                                          Population == "Ouahigouya" ~ "Africa",
                                                                          Population == "Santarem" ~ "Out_of_Africa",
                                                                          Population == "Bangkok" ~ "Out_of_Africa",
                                                                          Population == "Bebedouro" ~ "Out_of_Africa",
                                                                          Population == "Tafuna_Village" ~ "Out_of_Africa",
                                                                          Population == "Zac_Panda" ~ "Out_of_Africa",
                                                                          Population == "Jeddah" ~ "Out_of_Africa",
                                                                          Population == "Samut_Sakhon" ~ "Out_of_Africa",
                                                                          Population == "Tapachula" ~ "Out_of_Africa",
                                                                          TRUE  ~ "warning")
                                          ) %>% as.data.frame();



# merge both datasets
Populations_diversity_metrics_stats <- rbind(nucdiv_stats.Population_groups,tajsD_stats.Population_groups);


# output file for ST-5: stats by population
write.table(Populations_diversity_metrics_stats, file = "aedes_aegypti.downsampled_populations.angsd_pi_tajD_scores_byPop.txt", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE);








# ------------------------------------------------------------------------------------------------- #
#                         *** Reproduce stats by groups: Table 1 ***                                #
# ------------------------------------------------------------------------------------------------- #

# Change "headers" name to "capital letters" to avoid conflict between functions and column names while making summary stats
names(nucdiv_stats.Population_groups) <- toupper(names(nucdiv_stats.Population_groups));
names(tajsD_stats.Population_groups)  <- toupper(names(tajsD_stats.Population_groups));


# get stats by groups

# nucleotide diversity (π)
# -------------------------------------------------------------------------------------------------
nucdiv_stats.PopGroups <- nucdiv_stats.Population_groups %>% 
                                          group_by(POPGROUP) %>% 
                                                   summarise(df = length(MEAN)-1,
                                                               min  = min(MEAN),
                                                               max  = max(MEAN),
                                                               mean = mean(MEAN),
                                                               sd   = sd(MEAN),
                                                               se   = sd(MEAN)/sqrt(length(MEAN)),
                                                               t.score = qt(p=alpha/2, df=df,lower.tail=F),
                                                               margin.error = t.score * se,
                                                               ci.lower = mean - margin.error,
                                                               ci.upper = mean + margin.error
                                                   ) %>% as.data.frame()




# Tajima's D
# -------------------------------------------------------------------------------------------------
tajsD_stats.PopGroups <- tajsD_stats.Population_groups %>% 
                                          group_by(POPGROUP) %>% 
                                                   summarise(df = length(MEAN)-1,
                                                               min  = min(MEAN),
                                                               max  = max(MEAN),
                                                               mean = mean(MEAN),
                                                               sd   = sd(MEAN),
                                                               se   = sd(MEAN)/sqrt(length(MEAN)),
                                                               t.score = qt(p=alpha/2, df=df,lower.tail=F),
                                                               margin.error = t.score * se,
                                                               ci.lower = mean - margin.error,
                                                               ci.upper = mean + margin.error
                                                   ) %>% as.data.frame()





# -------------------------------------------------------------------------------------------------
# Collect the data:
# Table 1 | Measures of genetic diversity for the sampled Ae. aegypti mosquitoes: nucleotide diverity (pi) and Tajima's D for Africa, Out of Africa, and Aaa-like populations. 
# -------------------------------------------------------------------------------------------------

# *** Global: nucleotide diversity ***
nucdiv_stats.PopGroups %>% dplyr::select(POPGROUP,mean);
#        POPGROUP        mean
#          Africa  0.03697758
#   Out_of_Africa  0.01436593


# *** Global: tajima's D ***
tajsD_stats.PopGroups %>% dplyr::select(POPGROUP,mean);
#        POPGROUP        mean
#          Africa  -0.8514954
#   Out_of_Africa   0.2778591


# *** nucleotide diversity + tajima's D Aaa-like pops ***
Populations_diversity_metrics_stats %>% dplyr::select(POPGROUP,Population,dataset,mean) %>% dplyr::filter(Population=="Thies" | Population=="Ngoye" | Population=="Rabai_dom")
#     POPGROUP Population dataset        mean
#       Africa      Thies  NucDiv  0.02248769
#       Africa      Ngoye  NucDiv  0.02157758
#       Africa  Rabai_dom  NucDiv  0.01616387
#       Africa      Thies    TajD -0.68677633
#       Africa      Ngoye    TajD -0.76077233
#       Africa  Rabai_dom    TajD -0.47383967



# Well, that's it! Bye!
```
