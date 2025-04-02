#!/usr/bin/env Rscript

# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada Chavez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #



# /// START ///
#########################################################################################################################################
Sys.time()

# set library source
.libPaths("/home/nabor/R/R-3.6.2/lib/R/library/");
rm(list=ls())


library(pcadapt); require(qvalue);  library(tibble);



# set directory
setwd("/directory/");


# load snps (plink format)
path_to_file <- "mySNPs.bed";

filename <- read.pcadapt(path_to_file, type = "bed");

# 'mahalanobis' approach:
system.time( pcadapt_outliers <- pcadapt(input = filename, K = 6, method = "mahalanobis", LD.clumping = list(size = 200, thr = 0.1), min.maf = 0.01) );

# transform pvalues to qvalues:
qval <- qvalue(pcadapt_outliers$pvalues)$qvalues;

# 1% of the dataset will be false positive
alpha <- 0.01;

outliers.qvalues <- which(qval < alpha);



# get SNP outliers positions (plink format)
path_to_BIM_file <- "mySNPs.bim";
mapped_positions <- bigreadr::fread2(path_to_BIM_file);

# Get information of outliers associated to their VCF positions and other basic info"
outliers_positions <- data.frame(Chrm      = mapped_positions$V1[outliers.qvalues], 
                                 Pos       = mapped_positions$V4[outliers.qvalues], 
                                 snpID     = mapped_positions$V2[outliers.qvalues], 
                                 REF       = mapped_positions$V6[outliers.qvalues],  ## this is REF = A2 (plink format)
                                 ALT       = mapped_positions$V5[outliers.qvalues],  ## this is ALT = A1 (plink format)
                                 vcfRowPos = outliers.qvalues);


# SNPs associated to each Principal Component: 
outliers_byPC <- get.pc(pcadapt_outliers, outliers.qvalues);

# merge full info
outliers_info <- cbind(outliers_positions, outliers_byPC);


# counter of positions
outliers_main_info <- tibble::rownames_to_column(outliers_info, "Nr");


# outfile
write.table(outliers_main_info, "outliers_outfile.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE );



# /// END ///
#########################################################################################################################################

