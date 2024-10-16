## Detection of selection using the "McDonald and Kreitman test (MKT)"

1. We first performed a search of orthologs genes between Aedes aegypti and Aedes albopictus across the whole set of 14,677 genes in the reference genome AaegL5 as described in Methods.
2. Next, for each single individual sample in an African and Out of Africa population, we created a whole genome sequences using the reference genome AaegL5, and reference nucleotides were replaced by the corresponding SNPs nucleotide change using a .
3. After downsampling each population (as described in Methods), we extracted all orthologs genes in each single individual genome and population.
4. Codon alignments were created and refined using 'macse' v2.07, and verified with 'pal2nal' program.
5. We used the 'iMKT' R package and scripts provided by the [Murga-Moreno etal 2019: iMTK group paper](https://academic.oup.com/nar/article/47/W1/W283/5488529?login=false) to detect signals of Selection on each single gene alignment in each population of Africa and Out of Africa.

```bash
## /// Step 1 ///

# Blast + clustering similarities:
proteinortho_grab_proteins.pl -tofiles protein_ortho_aaegL5_vs_alboFs.proteinortho.tsv  'Aedes-aegypti-LVP_AGWG_AaegL5_2.longest_isoforms.faa'  'VectorBase-55_AalbopictusFoshanFPA.longest_isoforms.faa'  'VectorBase-61_AalbopictusFoshan.longest_isoforms.faa'  -p=blastp+  -cpus=60  -sim=1  -18 singles  -xml  -identity=0.25  -coverage=50  -evalue=0.00001

# Orthology
proteinortho_clustering  protein_ortho_aaegL5_vs_alboFs.blast-graph  >  protein_ortho_aaegL5_vs_alboFs.proteinortho-graph.main_output.txt


## /// Step 2 to 5 ///
```bash
$CDS_POPS=/home/nlozada/aaegypti/all_downsampled_pops/alignments/*_CDS

for i in CDS_POP;
   pop_fasta_genome = ${i}.fasta;
   pop_VCF_name = ${i}.vcf.gz;
   bcftools consensus  --fasta-ref Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta  --output $pop_fasta_genome  $pop_VCF_name;
   agat_sp_extract_sequences.pl -g Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.GFF3 -f $pop_fasta_genome -type cds;
done


## /// Step 4 ///

# Downsampled populations
# list of populations whole genome coding protein genes: main DIRs 
ALL_POPS_DIR_PATH_NAMES=full_path_whole_proteome_all_poppulations_DIRnames.txt
ALL_POPS_PROTEOME_GENEID=full_path_whole_proteome_population_GENEID.txt;

for POPDIR in  ALL_POPS_DIR_PATH_NAMES;
   cd $POPDIR;
   for GENEID in ALL_POPS_PROTEOME_GENEID;
      java -jar macse.jar $GENEID   
      java -jar macse_v2.07.jar  -prog alignSequences  -seq ${GENEID}.cds.ortho.snps.fasta  -max_refine_iter 4  -out_AA ${GENEID}.cds.ortho.snps.AA.aln.fasta  -out_NT ${GENEID}.cds.ortho.snps.NT.aln.fasta;
      java -jar macse_v2.07.jar  -prog exportAlignment  -align ${GENEID}.cds.ortho.snps.NT.aln.fasta   -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS ---  -out_AA ${GENEID}.cds.ortho.snps.AA.aln_noFS.fasta  -out_NT ${GENEID}.cds.ortho.snps.NT.aln_noFS.fasta;
      # sort alignment as requested in "Murga-Moreno etal 2019: iMTK group paper": 1) reference gene ID AaegL5, 2) all genes (same geneID as 1) with the polymorphic mutations (SNPs), 3) ortholog gene
      # next, validate codon alignment
      pal2nal.pl  ${GENEID}.cds.ortho.snps.AA.aln_noFS.sorted.fasta  ${GENEID}.cds.ortho.snps.NT.aln_noFS.sorted.fasta  -nomismatch -nogap -output fasta > ${GENEID}.cds.ortho.snps.AA.aln_noFS.sorted.codon.fasta;
      done
done

## /// Step 5 ///
# calculate DAF and DIV:

ALL_POPS_CODON_ALNS_GENEIDS=FULLPATH_POP_CODON_ALNS_GENEID.txt;

# changes columns from the original DAF and DIV files separately;
parser_command_string_DAF = "'{print \$1\"\\t\"\$3\"\\t\"\$2}'";
parser_command_string_DIV = "'{print \$3\"\\t\"\$2\"\\t\"\$4\"\\t\"\$1}'";


foreach CODON in ALL_POPS_CODON_ALNS_GENEIDS;
   CODON_name = FILE2=$(basename ${FILE1} .fasta)
   do python2.7 sfsFromFasta_v2.py --multiFasta $CODON  --daf {CODON}.daf  --div {CODON}.div  --codonTable standard
   awk '{print $1,"\t",$3,"\t",$2}' {CODON}.daf > {CODON}.parsed.daf;
   awk '{print $1,"\t",$3,"\t",$2}' {CODON}.div > {CODON}.parsed.div;
done
```

```R
# for each gene ID in each single populoatioon, calculate the MKT and DoS:

# read files with DAF names, and also read DIV ones.
for(i in 1:length(dafs)){
    daf_file <- paste0(dafs[i]);
    div_file <- gsub("fas2daf.parsed.daf","fas2div.parsed.div",daf_file);
    daf      <- try(read.delim(daf_file, header = T, stringsAsFactors = F));
    div      <- try(read.delim(div_file, header = T, stringsAsFactors = F));

    mkt_list[[i]]  <- try(standardMKT(daf, div));
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

myGENE_ID <- gsub('.*\\/', '', dafs);

# MK-test
MKT_STD <- data.frame(geneID=myGENE_ID %>% gsub(pattern=".cds.downsampled.ortho.snps.codon.NT.aln.noFS.sorted.fas2daf.parsed.daf", replace=""),
                          MKT_STD_Fishers_pvalue=mkt_list %>% lapply("[",2) %>% unlist(),
                          MKT_std_alpha=mkt_list %>% lapply("[",1) %>% unlist(),
                          stringsAsFactors = FALSE) #%>% subset(!is.na(MKT_STD_pvalue))
  

MKT_STD$MKT_STD_alpha_clean <- gsub("\\n", "#", MKT_STD$MKT_std_alpha)

# MK-test standard
MKT_STD$MKT_STD_Fisher_pval.adj <- p.adjust(MKT_STD$MKT_STD_Fishers_pvalue, method = "BH");



# DGRP
MKT_DGRP <- data.frame(geneID2=myGENE_ID %>% gsub(pattern=".cds.downsampled.ortho.snps.codon.NT.aln.noFS.sorted.fas2daf.parsed.daf",replace=""),
                          MKT_STD_alpha=DGRP_list %>% lapply(getAlpha) %>% unlist(),
                          MKT_div_Ka=DGRP_list %>% lapply(getDivMetrics_Ka) %>% unlist(),
                          MKT_div_Ks=DGRP_list %>% lapply(getDivMetrics_Ks) %>% unlist(),
                          MKT_div_omega=DGRP_list %>% lapply(getDivMetrics_OMEGA) %>% unlist(),
                          MKT_MKtbl_pS_counts=DGRP_list %>% lapply(getMKcontigencyTbl_NEUTRAL_Polymorphism) %>% unlist(),
                          MKT_MKtbl_dS_counts=DGRP_list %>% lapply(getMKcontigencyTbl_NEUTRAL_Divergence) %>% unlist(),
                          MKT_MKtbl_pN_counts=DGRP_list %>% lapply(getMKcontigencyTbl_SELECTED_Polymorphism) %>% unlist(),
                          MKT_MKtbl_dN_counts=DGRP_list %>% lapply(getMKcontigencyTbl_SELECTED_Divergence) %>% unlist(),
                          MKT_DGRP_fraction_d_cutoff0=DGRP_list %>% lapply(getFractions_0d) %>% unlist(),
                          MKT_DGRP_fraction_f_cutoff0=DGRP_list %>% lapply(getFractions_0f) %>% unlist(),
                          MKT_DGRP_fraction_b_cutoff0=DGRP_list %>% lapply(getFractions_0b) %>% unlist(),
                          MKT_DGRP_fraction_d_cutoff005=DGRP_list %>% lapply(getFractions_005d) %>% unlist(),
                          MKT_DGRP_fraction_f_cutoff005=DGRP_list %>% lapply(getFractions_005f) %>% unlist(),
                          MKT_DGRP_fraction_b_cutoff005=DGRP_list %>% lapply(getFractions_005b) %>% unlist(),
                          MKT_DGRP_fraction_d_cutoff01=DGRP_list %>% lapply(getFractions_01d) %>% unlist(),
                          MKT_DGRP_fraction_f_cutoff01=DGRP_list %>% lapply(getFractions_01f) %>% unlist(),
                          MKT_DGRP_fraction_b_cutoff01=DGRP_list %>% lapply(getFractions_01b) %>% unlist(),
                          stringsAsFactors = FALSE)



# merge all MK test tables
main_MKT_results <- cbind(MKT_STD,MKT_DGRP);



# Calculate DoS
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

