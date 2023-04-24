#!/usr/bin/env Rscript

# -------------------------------------------- #
# Questions:                                   #                        
# Alejandro Nabor Lozada Chavez                #        
# nabor.lozada@gmail.com                       #           
# alejandro.chavez@unipv.it                    #        
# -------------------------------------------- #

# /// START ///
#########################################################################################################################################
Sys.time()


# set library source:
.libPaths("/home/jaketu/R/R-3.6.2/lib/R/library/")
rm(list=ls())



# Required libraries
library(GenomicFeatures);   library(Rsamtools);   library(VariantAnnotation);   library(tidyverse);


# Set directory
setwd("/directory/");



# // Annotations: Refenrece geneme (fa), gene annotations (txdb), SNPs (vcf) //
# reference fasta genome
fa=open(FaFile("Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa"));
# reference gff genome
txdb=makeTxDbFromGFF(file="Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.gff3", format="gff3");
# target VCF
vcf=readVcf("mySNPs.vcf.gz");



# get effects in Protein Coding Genes
effects = predictCoding(vcf, txdb, fa);



# list of genes in the genome/gff
AaegL5_genes = unique(effects$GENEID);

# define variable
main_results_table=data.frame(GeneID=AaegL5_genes,L0=0,L2=0,L4=0,A0=0,A2=0,A4=0,B0=0,B2=0,B4=0);

# Effects by gene ID
geneEffects=lapply(AaegL5_genes, function(x,effects) effects[effects$GENEID==x], effects=effects);



# Get the reference codons (REFCODONS) for each gene, and variant codons (VARCODONS)
rC = lapply(geneEffects, function(x) unname(as.vector(x$REFCODON, mode="character")));
vC = lapply(geneEffects, function(x) unname(as.vector(x$VARCODON, mode="character")));


# only codons (3nts)
for (i in seq_along(rC)) {
    removeThis = unique(c(which(nchar(rC[[i]])>3), which(nchar(vC[[i]])>3)))
    if (length(to.remove)>0) {
        rC[[i]]=rC[[i]][-removeThis]
        vC[[i]]=vC[[i]][-removeThis]
    }
}





# Function to classify each site in the gene as nondegenerate, 2-fold degen, or 4-fold degen
# ----------------------------------------------------------------------------------------------------------------------------------------
table_degenerates_sites <- function(codons_list) {
                                    degeneracy=list(T=list(T=list(T=c(3,3,2), C=c(3,3,2), A=c(2,3,2), G=c(2,3,2)),
                                                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                                                           A=list(T=c(3,3,2), C=c(3,3,2), A=c(3,2,2), G=c(3,3,2)),
                                                           G=list(T=c(3,3,2), C=c(3,3,2), A=c(3,2,3), G=c(3,3,3))),
                                                    C=list(T=list(T=c(2,3,0), C=c(3,3,0), A=c(2,3,0), G=c(2,3,0)),
                                                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                                                           A=list(T=c(3,2,2), C=c(3,2,2), A=c(3,2,2), G=c(3,2,2)),
                                                           G=list(T=c(3,3,0), C=c(3,3,0), A=c(2,3,0), G=c(2,3,0))),
                                                    A=list(T=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                                                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                                                           A=list(T=c(3,2,2), C=c(3,2,2), A=c(3,2,2), G=c(3,2,2)),
                                                           G=list(T=c(3,3,2), C=c(3,3,2), A=c(2,3,2), G=c(2,3,2))),
                                                    G=list(T=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                                                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                                                           A=list(T=c(3,2,2), C=c(3,2,2), A=c(3,2,2), G=c(3,2,2)),
                                                           G=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0))))
                                    d=unlist(lapply(codons_list,function(x) eval(parse(text=paste0(c("degeneracy",unlist(strsplit(x, split="")),use.names=F), collapse="$")))))
                                    d=as.numeric(gsub("3", "0", gsub(0,4,d)))
                                    return(d)
                            }



# Convert each set of rC to a vector of 0s, 2s, and 4s, to indicate level of degeneracy for every site
degC=lapply(rC, table_degenerates_sites)



# Calculate L0, L2, L4
main_results_table$L0 = unlist(lapply(degC, function(x) (length(which(x==0)))/(length(x))));
main_results_table$L2 = unlist(lapply(degC, function(x) (length(which(x==2)))/(length(x))));
main_results_table$L4 = unlist(lapply(degC, function(x) (length(which(x==4)))/(length(x))));



# Transitins(Ts) or Transversions(Tv)
TsTv <- function(rC, vC) {
                    ref=paste0(rC, collapse="")
                    var=paste0(vC, collapse="")
                    ref=unlist(strsplit(ref, split=""), use.names=F)
                    var=unlist(strsplit(var, split=""), use.names=F)
                    mut = which(ref != var)
                    temp=ref[mut]
                    temp[which(ref[mut]=="G")] = "A"
                    temp[which(ref[mut]=="A")] = "G"
                    temp[which(ref[mut]=="C")] = "T"
                    temp[which(ref[mut]=="T")] = "C"
                    transV = which(temp != var[mut])
                    res = rep("N", length(ref))
                    res[mut[transV]] = "Tv"
                    res[mut[-transV]] = "Ts"
                    return(res)
        }


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


# Method's Li (1993): improved Kimura's 2-parameter model to calculate pN and pS
main_results_table$pS = main_results_table$B4 + (((main_results_table$L2 * main_results_table$A2) + (main_results_table$L4*main_results_table$A4))/(main_results_table$L2 + main_results_table$L4))
main_results_table$pN = main_results_table$A0 + (((main_results_table$L0 * main_results_table$B0) + (main_results_table$L2*main_results_table$B2))/(main_results_table$L0 + main_results_table$L2))

# reference: J. Mol. Evol. 1993: 36,96-99. 
# Equations 8 and 9, respectively.


# And finally, the pNpS ratio:
main_results_table$pNpS = main_results_table$pN/main_results_table$pS




neutrality = 1.0;

my_results_pnps <- main_results_table %>% mutate(selection = case_when(pNpS > neutrality  ~ "positive",
                                                                       pNpS < neutrality  ~ "negative",
                                                                       pNpS == neutrality ~ "neutral",
                                                                       pNpS == "Inf"      ~ "WARNING",
                                                                       TRUE  ~ "background"));



# Save main table with results
write.table(my_results_selection, "outfile.stats.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE );

# /// END ///
#########################################################################################################################################

