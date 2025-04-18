## Detection of selection using the "McDonald and Kreitman test (MKT)"

1. We first performed a search of orthologs genes between Aedes aegypti and Aedes albopictus across the whole set of 14,677 genes in the reference genome AaegL5 as described in Methods.
2. Next, for each single individual sample in an African and Out of Africa population, we created a whole genome sequences using the reference genome AaegL5, and reference nucleotides were replaced by the corresponding SNPs nucleotide change using a .
3. After downsampling each population (as described in Methods), we extracted all orthologs genes in each single individual genome and population.
4. Codon alignments were created and refined using 'macse' v2.07, and verified with 'pal2nal' program.
5. We used the 'iMKT' R package and scripts provided by the [Murga-Moreno etal 2019: iMTK group paper](https://academic.oup.com/nar/article/47/W1/W283/5488529?login=false) to detect signals of Selection on each single gene alignment in each population of Africa and Out of Africa.
6. Calculate the Direction of Selection (DoS) for each gene ID of each population.

```bash

## /// Step 1 ///
## --------------------------------------------------
cd /home/nlozada/aaegypti/detection_selection/;
mkdir ortholgs.aaeg_vs_aalbo; 
cd ortholgs.aaeg_vs_aalbo; 

# Blast + clustering similarities:
proteinortho_grab_proteins.pl -tofiles protein_ortho_aaegL5_vs_alboFs.proteinortho.tsv  'Aedes-aegypti-LVP_AGWG_AaegL5_2.longest_isoforms.faa'  'VectorBase-55_AalbopictusFoshanFPA.longest_isoforms.faa'  'VectorBase-61_AalbopictusFoshan.longest_isoforms.faa'  -p=blastp+  -cpus=60  -sim=1  -18 singles  -xml  -identity=0.25  -coverage=50  -evalue=0.00001

# Orthology
# process network *.blast-graph file
proteinortho_clustering  protein_ortho_aaegL5_vs_alboFs.blast-graph  >  protein_ortho_aaegL5_vs_alboFs.proteinortho-graph.main_output.txt


## /// Step 2 to 5 ///

# population directories and single sample 'reconstructed pseudogenomes' with SNPs.
$CDS_POPS=/home/nlozada/aaegypti/all_downsampled_pops/alignments/*_CDS

for i in CDS_POP; do
    cd $i;
    for j in *.fasta; do
        pop_fasta_genome = ${j}.fasta;
        pop_VCF_name = ${j}.vcf.gz;
        bcftools consensus  --fasta-ref Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta  --output $pop_fasta_genome  $pop_VCF_name;
        agat_sp_extract_sequences.pl -g Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.GFF3 -f $pop_fasta_genome -type cds;
        cd ..;
    done
done



## /// Step 4 ///
## --------------------------------------------------
# Downsampled populations
# list of populations whole genome coding protein genes: main DIRs 
ALL_POPS_DIR_PATH_NAMES=full_path_whole_proteome_all_poppulations_DIRnames.txt
ALL_POPS_PROTEOME_GENEID=full_path_whole_proteome_population_GENEID.txt;

for POPDIR in  ALL_POPS_DIR_PATH_NAMES;
    cd $POPDIR;
    for GENEID in ALL_POPS_PROTEOME_GENEID;
        java -jar macse_v2.07.jar  -prog alignSequences  -seq ${GENEID}.cds.ortho.snps.fasta  -max_refine_iter 4  -out_AA ${GENEID}.cds.ortho.snps.AA.aln.fasta  -out_NT ${GENEID}.cds.ortho.snps.NT.aln.fasta;
        java -jar macse_v2.07.jar  -prog exportAlignment  -align ${GENEID}.cds.ortho.snps.NT.aln.fasta   -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS ---  -out_AA ${GENEID}.cds.ortho.snps.AA.aln_noFS.fasta  -out_NT ${GENEID}.cds.ortho.snps.NT.aln_noFS.fasta;
        # sort alignment as requested in "Murga-Moreno etal 2019: iMTK group paper": 1) reference gene ID AaegL5, 2) all genes (same geneID as 1) with the polymorphic mutations (SNPs), 3) ortholog gene
        # next, validate codon alignment
        pal2nal.pl  ${GENEID}.cds.ortho.snps.AA.aln_noFS.sorted.fasta  ${GENEID}.cds.ortho.snps.NT.aln_noFS.sorted.fasta  -nomismatch -nogap -output fasta > ${GENEID}.cds.ortho.snps.AA.aln_noFS.sorted.codon.fasta;
        done
    cd ..;
done



## /// Step 5 ///
## --------------------------------------------------

# calculate DAF and DIV:

ALL_POPS_CODON_ALNS_GENEIDS=FULLPATH_POP_CODON_ALNS_GENEID.txt;

# then, changes columns order from the original DAF and DIV files separately;

fore CODON in ALL_POPS_CODON_ALNS_GENEIDS;
   CODON_name = $(basename ${CODON} .fasta)
   do python2.7 sfsFromFasta.py --fasta $CODON  --daf {CODON}.daf  --div {CODON}.div  --codonTable standard
   awk '{print $1,"\t",$3,"\t",$2}' {CODON}.daf > {CODON}.parsed.daf;
   awk '{print $1,"\t",$3,"\t",$2}' {CODON}.div > {CODON}.parsed.div;
done
```

## Calculation of MKT and DoS for each gene and population.

This sample custom R script is used for all genes to be analyzed in a single population. A list of genes most be provided so the R script can read it and parse the names in it to read each single pair of infiles: `daf` and `div`. Here we used the `daf` files as reference. The extension names of the gene files can be removed and/or replace in the code.


```R
# for each gene ID in each single populoatioon, calculate the MKT and DoS:


suppressMessages(library(iMKT));
suppressMessages(library(dplyr));
suppressMessages(library(tidyr));

# File names of DAF calculations
DAF_FILES <- readLines("/home/nlozada/aedes_aegypti/NEE_paper/results/make_MKtest/daf_n_div_files/country.popname1.downsampled.daf_list_files.txt")
# pop 1 out of 40

# read files with DAF names, parsed names to read also DIV files. Make MKT:
for(i in 1:length(DAF_FILES)){
    daf_file <- paste0(DAF_FILES[i]);
    div_file <- gsub("fas2daf.parsed.daf","fas2div.parsed.div",daf_file);
    daf      <- try(read.delim(daf_file, header = T, stringsAsFactors = F));
    div      <- try(read.delim(div_file, header = T, stringsAsFactors = F));

    mkt_list[[i]]   <- try(standardMKT(daf, div));
    DGRP_list[[i]]  <- try(DGRP(daf, div));
}



# ------------------------------------------------------------------------------
# functions to extract values:
# ------------------------------------------------------------------------------

# pvalue
getPval <- function(x){
  if(!inherits(x,"try-error")){
    x$Results$`Fishers exact test P-value`[2] %>% return()
  } else {
    return(NA)
  }
}

# get alpha from Fisher test
getAlpha <- function(x){
  if(!inherits(x,"try-error")){
    x$Results$alpha.symbol[2] %>% return()
  } else {
    return(NA)
  }
}



# divergence metrics
# -------------------------------------------------------------------------------
# nonsyn
getDivMetrics_Ka <- function(x){
  if(!inherits(x,"try-error")){ x$`Divergence metrics`$`Global metrics`[1] %>% return() } else { return(NA) }
}

# synon
getDivMetrics_Ks <- function(x){
  if(!inherits(x,"try-error")){ x$`Divergence metrics`$`Global metrics`[2] %>% return() } else { return(NA) }
}


# MK test table
# -------------------------------------------------------------------------------

getMKcontigencyTbl_NEUTRAL_Polymorphism <- function(x){
  if(!inherits(x,"try-error")){ x$`MKT tables`$`MKT standard table`[1,1] %>% return() } else { return(NA) }
}

getMKcontigencyTbl_NEUTRAL_Divergence <- function(x){
  if(!inherits(x,"try-error")){ x$`MKT tables`$`MKT standard table`[1,2] %>% return() } else { return(NA) }
}

getMKcontigencyTbl_SELECTED_Polymorphism <- function(x){
  if(!inherits(x,"try-error")){ x$`MKT tables`$`MKT standard table`[2,1] %>% return() } else { return(NA) }
}

getMKcontigencyTbl_SELECTED_Divergence <- function(x){
  if(!inherits(x,"try-error")){ x$`MKT tables`$`MKT standard table`[2,2] %>% return() } else { return(NA) }
}




# -------------------------------------------------------------------------------
# MKT: Finally get main outputs
# -------------------------------------------------------------------------------

myGENE_ID <- gsub('.*\\/', '', DAF_FILES);


MKT_STD <- data.frame(geneID=myGENE_ID %>% gsub(pattern=".cds.downsampled.ortho.snps.codon.NT.aln.noFS.sorted.fas2daf.parsed.daf", replace=""),
                          MKT_STD_Fishers_pvalue=mkt_list %>% lapply("[",2) %>% unlist(),
                          MKT_std_alpha=mkt_list %>% lapply("[",1) %>% unlist(),
                          stringsAsFactors = FALSE) #%>% subset(!is.na(MKT_STD_pvalue))


# extract alpha scores from MK-test standard
MKT_STD$MKT_STD_alpha_clean <- gsub("\\n", "#", MKT_STD$MKT_std_alpha)

# extract adjusted pval MK-test standard
MKT_STD$MKT_STD_Fisher_pval.adj <- p.adjust(MKT_STD$MKT_STD_Fishers_pvalue, method = "BH");



# DGRP (same as MKT standard, but with counts for neutral and deleterious mutations proportion on the protein coding gene)
MKT_DGRP <- data.frame(geneID2=myGENE_ID %>% gsub(pattern=".cds.downsampled.ortho.snps.codon.NT.aln.noFS.sorted.fas2daf.parsed.daf",replace=""),
                          MKT_STD_alpha=DGRP_list %>% lapply(getAlpha) %>% unlist(),
                          MKT_div_Ka=DGRP_list %>% lapply(getDivMetrics_Ka) %>% unlist(),
                          MKT_div_Ks=DGRP_list %>% lapply(getDivMetrics_Ks) %>% unlist(),
                          MKT_MKtbl_pS_counts=DGRP_list %>% lapply(getMKcontigencyTbl_NEUTRAL_Polymorphism) %>% unlist(),
                          MKT_MKtbl_dS_counts=DGRP_list %>% lapply(getMKcontigencyTbl_NEUTRAL_Divergence) %>% unlist(),
                          MKT_MKtbl_pN_counts=DGRP_list %>% lapply(getMKcontigencyTbl_SELECTED_Polymorphism) %>% unlist(),
                          MKT_MKtbl_dN_counts=DGRP_list %>% lapply(getMKcontigencyTbl_SELECTED_Divergence) %>% unlist(),
                          stringsAsFactors = FALSE)


# merge all MK test tables
main_MKT_results <- cbind(MKT_STD,MKT_DGRP);



# /// Step 6: Calculate DoS ///
# -------------------------------------------------------------------------------

# Direction of selection:
# DoS = Dn/( Dn + Ds ) - Pn/( Pn + Ps )
# 
# Pn: nonsilent polymorphisms
# Ps: silent polymorphisms
# Dn: nonsilent and 
# Ds: silent substitutions

# DOS dnds component: D(n)/( D(n) + D(s) ) 
main_MKT_results$DOS_dNdS_STAT <- main_MKT_results$MKT_MKtbl_dN_counts / (main_MKT_results$MKT_MKtbl_dN_counts + main_MKT_results$MKT_MKtbl_dS_counts);

# DOS pnps component: P(n)/( P(n) + P(s) )
main_MKT_results$DOS_pNpS_STAT <- main_MKT_results$MKT_MKtbl_pN_counts / (main_MKT_results$MKT_MKtbl_pN_counts + main_MKT_results$MKT_MKtbl_pS_counts);

main_MKT_results$DOS_STAT <- main_MKT_results$DOS_dNdS_STAT - main_MKT_results$DOS_pNpS_STAT;


main_MKT_results$total_Polymorphic_sites <- main_MKT_results$MKT_MKtbl_pS_counts + main_MKT_results$MKT_MKtbl_pN_counts;
main_MKT_results$total_Divergent_sites   <- main_MKT_results$MKT_MKtbl_dS_counts + main_MKT_results$MKT_MKtbl_dN_counts;


write.table(main_MKT_results, "OUTFILE_NAME.MKT_n_DoS.TXT", quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE )



```

