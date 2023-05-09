#!/usr/bin/env Rscript

# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada Chavez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #




# Li's Method (1993): improved Kimura's 2-parameter model
# reference: J. Mol. Evol. 1993: 36,96-99. 




# /// START ///
#########################################################################################################################################
Sys.time()


# set library source:
.libPaths("/home/nabor/R/R-3.6.2/lib/R/library/")
rm(list=ls())



# Required libraries
library(GenomicFeatures);   library(Rsamtools);   library(VariantAnnotation);   library(tidyverse);


# Set directory
setwd("/directory/");



# // Annotations //
# reference fasta genome
fa=open(FaFile("Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa"));
# reference gff genome
txdb=makeTxDbFromGFF(file="Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.gff3", format="gff3");
# target VCF
vcf=readVcf("mySNPs.vcf.gz");



# get effects
effects = predictCoding(vcf, txdb, fa);

# list of genes
AaegL5_genes = unique(effects$GENEID);

# define variable
main_results_table=data.frame(GeneID=AaegL5_genes,L0=0,L2=0,L4=0,A0=0,A2=0,A4=0,B0=0,B2=0,B4=0);

# Effects
geneEffects=lapply(AaegL5_genes, function(x,effects) effects[effects$GENEID==x], effects=effects);

# Get the reference codons (REFCODONS) and variant codons (VARCODONS)
rC = lapply(geneEffects, function(x) unname(as.vector(x$REFCODON, mode="character")));
vC = lapply(geneEffects, function(x) unname(as.vector(x$VARCODON, mode="character")));



# Level of degeneracy: rC to a vector of 0s, 2s, and 4s
degC=lapply(rC, table_degenerates_sites)


# Calculate L0, L2, L4
main_results_table$L0 = unlist(lapply(degC, function(x) (length(which(x==0)))/(length(x))));
main_results_table$L2 = unlist(lapply(degC, function(x) (length(which(x==2)))/(length(x))));
main_results_table$L4 = unlist(lapply(degC, function(x) (length(which(x==4)))/(length(x))));


# find Ts, Tv
TsTv_byGene = vector("list", length(rC))


for (i in seq_along(rC)) {
    TsTv_byGene[[i]]=TsTv(rC[[i]], vC[[i]])
}


# Counts for A0, A2, A4, B0, B2,and B4 for each gene
for (i in seq_along(TsTv_byGene)) {
    tmp=TsTv_byGene[[i]]
    x=degC[[i]]
    main_results_table$A0[i]=length(intersect(which(tmp=="Ts"), which(x==0)))
    main_results_table$A2[i]=length(intersect(which(tmp=="Ts"), which(x==2)))
    main_results_table$A4[i]=length(intersect(which(tmp=="Ts"), which(x==4)))
    main_results_table$B0[i]=length(intersect(which(tmp=="Tv"), which(x==0)))
    main_results_table$B2[i]=length(intersect(which(tmp=="Tv"), which(x==2)))
    main_results_table$B4[i]=length(intersect(which(tmp=="Tv"), which(x==4)))
}


# rates
main_results_table$pS = main_results_table$B4 + (((main_results_table$L2 * main_results_table$A2) + (main_results_table$L4*main_results_table$A4))/(main_results_table$L2 + main_results_table$L4))
main_results_table$pN = main_results_table$A0 + (((main_results_table$L0 * main_results_table$B0) + (main_results_table$L2*main_results_table$B2))/(main_results_table$L0 + main_results_table$L2))


# And finally, the pNpS ratio:
main_results_table$pNpS = main_results_table$pN/main_results_table$pS


# Save main table with results
write.table(my_results_selection, "outfile.stats.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE );

# /// END ///
#########################################################################################################################################

