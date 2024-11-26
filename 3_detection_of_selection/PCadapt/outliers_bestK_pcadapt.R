#!/usr/bin/env Rscript


# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada Chavez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #

library(pcadapt);


# set directory
setwd("/directory/");


# load snps (plink format)
path_to_file <- "mySNPs.bed";

filename <- read.pcadapt(path_to_file, type = "bed");

# find best K with mahalanobis
system.time( clusters <- pcadapt(input = filename, K = 20, method = "mahalanobis", LD.clumping = list(size = 200, thr = 0.1), min.maf = 0.01) )


# screeplot
pdf("outfile_screeplot.pdf", height = 8.0, width = 10.0,useDingbats=FALSE)
   plot(clusters, option = "screeplot");
dev.off()

