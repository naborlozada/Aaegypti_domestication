#!/usr/bin/env Rscript


# Description
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# RUN R using version 4.3.1:
# /path/to/R_binaries/R-4.3.1/bin/R
#
#
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


rm(list=ls());

# R version 4.3.1:
.libPaths("/path/to/R/library/R-4.3.1/library/");


suppressPackageStartupMessages( library(tibble) );
suppressPackageStartupMessages( require(qvalue) );
suppressPackageStartupMessages( library(dplyr) );
suppressPackageStartupMessages( library(optparse) );
suppressPackageStartupMessages( library(qqman) );
suppressPackageStartupMessages( library(readr) );




# define options
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
option_list = list(
    make_option(c("-i", "--infile1"), type="character", default=NULL, help="ORIGINAL RAiSD LIST Report file (TXT file)", metavar="RAiSD_ReportList.pop_or_chrm_names.txt")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if ( is.null(opt$infile1) ){
  cat("\n\nDescription:
        This R script parses the output of RAiSD to detect `hard` selective sweeps based on outlier detection (top 0.05% & 1.0% tails) in a chromosome.\n\n")
  #print_help(opt_parser);

  stop(" ---> The only mandatory argument is `infile1`: RAiSD_ReportList.pop_or_chrm_names.txt
  \nUSAGE:
  \tRscript  selective_sweeps.raisd.LOCAL.get_tail_tops_1_n_0.05.R  -i infile1 \n\n", call.=FALSE)
}
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




cat("\n\n");
Sys.time();
cat("\n");

# options
setwd("/path/to/outputs/");

cat("WORKING DIR:");
getwd();
cat("\n");
cat(">>> infile: ", opt$infile1, "\n");









# functions to extract raisd output
# ---------------------------------------------------------------------------------------------------------
get_raisd_output <- function(Report, ...){
  filename = strsplit(Report, '/', fixed = TRUE)[[1]][12];
  chromTAG <- strsplit(filename, '.', fixed = TRUE)[[1]][20];
  chromTAG <- as.numeric(chromTAG);

  dataset = read.table(paste(Report[1]), header=F, skip=0)[,1:7];
  names(dataset) <- c("position","startW","endW","VAR","SFS","LD","mu_stat");
  dataset$chr    <- chromTAG;
  dataset_sorted <- dataset[c(8, 1, 2, 3, 4, 5, 6, 7)]
  return(dataset_sorted);
}












# ==================================================================================================================================================================================== #
#                                                                        *** START: PARSE RAISD INFILE ***                                                                             #
# ==================================================================================================================================================================================== #



INFILE = opt$infile1;
DIR    = gsub('/RAiSD_ReportList..txt','',INFILE);
myData = read_lines(INFILE, n_max = 1);

# merge dir+file
# full_path_Report = "/path/to/infile/with/RAiSD/results/RAiSD_Report.infle_name.chrm";
full_path_Report <-paste0(as.character(DIR),"/", as.character(myData))

cat("\n>>> target file: ", full_path_Report, "\n\n");

# V1 = Genomic location
# V2 = start (of the window)
# V3 = end (of the window)
# V4 = VAR
# V5 = SFS
# V6 = LD
# V7 = mu statistic

# > head( read.table(paste(full_path_Report[1]), header=F, skip=0)[,1:7], n=40 ); 
#       V1   V2    V3    V4    V5     V6     V7
# 1  16799   53 33545 4.461 1.233 0.2857  1.572
# 2  16855  127 33583 4.456 1.233 0.2885  1.585
# 3  17023  451 33595 4.415 1.233 0.6667  3.629
# 4  17071  511 33630 4.411 1.279 0.8247  4.652
# 5  17106  527 33685 4.417 1.324 1.0310  6.029




# get a core filename:
CORE_FILE_NAME <- sapply(strsplit(myData, ".", fixed = TRUE), function(x){paste(x[[2]], x[[3]], x[[16]], sep=".")});


# get data
myReport_full  <- get_raisd_output(full_path_Report);
head(myReport_full);

# cutoff = 99.95% ---> top 0.05%
THRESHOLD_05 <- 0.9995;
# cutoff = 99.0%   ---> top 1%
THRESHOLD_1 <- 0.990;

# 0.05%
thres_05 <- as.numeric(THRESHOLD_05);
topQ_05  <- thres_05*100;
threshold_05 <- quantile(x=myReport_full$mu_stat, probs = thres_05);

# 1.0%
thres_1 <- as.numeric(THRESHOLD_1);
topQ_1  <- thres_1*100;
threshold_1 <- quantile(x=myReport_full$mu_stat, probs = thres_1);

cat("\nFinished parsing RAISD outfile...\n\n");







# ==================================================================================================================================================================================== #
#                                                                 *** OUTLIERS EXTRACTION: 0.05% and 1.0% ***                                                                          #
# ==================================================================================================================================================================================== #




# get the top 0.05%
# -----------------------------------------------------------------------------
myReport_fltrd_A <- myReport_full %>% mutate(outlier_thr_05 = ifelse(mu_stat > threshold_05, "outlier", "background"));

myReport_fltrd_B <- myReport_fltrd_A %>% mutate(outlier_thr_1 = ifelse(mu_stat > threshold_1, "outlier", "background"));


# counts
cat("Total outliers 0.05%:");
myReport_fltrd_B %>% group_by(outlier_thr_05) %>% tally();
cat("\n\nTotal outliers 1.0%:");
myReport_fltrd_B %>% group_by(outlier_thr_1) %>% tally();

# get outliers of each cutoff
myReport_fltrd_outliers_thres_05 <- subset(myReport_fltrd_B, outlier_thr_05 == "outlier") %>% select(1:9);
myReport_fltrd_outliers_thres_1  <- subset(myReport_fltrd_B, outlier_thr_1 == "outlier") %>% select(1:8,10);


# make output 1
mytime <- format(Sys.time(), "v%d%b%Y");
OUTFILE_1 <- paste("selective_sweeps.raisd.",CORE_FILE_NAME,".full_dataset.outliers_n_background.", mytime,".txt", sep="") 
write.table(myReport_fltrd_B, OUTFILE_1, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE );


# make output 2
OUTFILE_2 <- paste("selective_sweeps.raisd.",CORE_FILE_NAME,".fltrd_dataset.outliers.thres_05.", mytime,".txt", sep="") 
write.table(myReport_fltrd_outliers_thres_05, OUTFILE_2, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE );


# make output 3
OUTFILE_3 <- paste("selective_sweeps.raisd.",CORE_FILE_NAME,".fltrd_dataset.outliers.thres_1.", mytime,".txt", sep="") 
write.table(myReport_fltrd_outliers_thres_1, OUTFILE_3, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE );

cat("\n\nFinished extracting outliers from RAISD `mu statistics` top 0.05% and 1.0% ...\n\n");








# ==================================================================================================================================================================================== #
#                                                                 *** OUTLIERS BED FILES:  0.05% and 1.0% ***                                                                          #
# ==================================================================================================================================================================================== #



# ---- cutoff 0.05% ---- 
myReport_fltrd_outliers_thres_05_BED <- myReport_fltrd_outliers_thres_05;
myReport_fltrd_outliers_thres_05_BED$end <- myReport_fltrd_outliers_thres_05$position;
colnames(myReport_fltrd_outliers_thres_05_BED)[2] <- "start";

myReport_fltrd_outliers_thres_05_BED <- myReport_fltrd_outliers_thres_05_BED[c(1, 2, 10, 3, 4, 5, 6, 7, 8, 9)]


OUTFILE_4 <- paste("/path/to/raisd/local_outputs/selective_sweeps.raisd.",CORE_FILE_NAME,".fltrd_dataset.outliers.thres_05.", mytime,".bed", sep="") 
write.table(myReport_fltrd_outliers_thres_05_BED, OUTFILE_4, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE );
OUTFILE_4



# ---- cutoff 1.0% ---- 
myReport_fltrd_outliers_thres_1_BED <- myReport_fltrd_outliers_thres_1;
myReport_fltrd_outliers_thres_1_BED$end <- myReport_fltrd_outliers_thres_1$position;
colnames(myReport_fltrd_outliers_thres_1_BED)[2] <- "start";

myReport_fltrd_outliers_thres_1_BED <- myReport_fltrd_outliers_thres_1_BED[c(1, 2, 10, 3, 4, 5, 6, 7, 8, 9)]


OUTFILE_5 <- paste("/path/to/raisd/local_outputs/selective_sweeps.raisd.",CORE_FILE_NAME,".fltrd_dataset.outliers.thres_1.", mytime,".bed", sep="") 
write.table(myReport_fltrd_outliers_thres_1_BED, OUTFILE_5, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE );
OUTFILE_5


cat("\nFinished making BED files of outliers with top 0.05% and 1.0%. Now, mapping to Aedes aegypti reference genome...\n\n");







# ==================================================================================================================================================================================== #
#                                                          *** MAP BED file to REFERENCE GENOME:  0.05% and 1.0% ***                                                                   #
# ==================================================================================================================================================================================== #




# mapping coordinates to genome
reference_genome <- "/path/to/reference_genome/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.bed";

# outfile
OUTFILE_6_thres_05_genome_v_outliers  <- paste("selective_sweeps.raisd.",CORE_FILE_NAME,".fltrd_dataset.outliers.thres_05.sweeps_to_refgenome.",mytime,".bed", sep="");
OUTFILE_7_thres_1_genome_v_outliers   <- paste("selective_sweeps.raisd.",CORE_FILE_NAME,".fltrd_dataset.outliers.thres_1.sweeps_to_refgenome.",mytime,".bed", sep="");


# ** top 0.05 % [REVIEWER 1] **
system( paste("time /path/to/software_scripts/bedtools intersect -a ", OUTFILE_4," -b", reference_genome," -wb | grep -v -w chromosome > ", OUTFILE_6_thres_05_genome_v_outliers) );


# ** top 1.0 % [ours] **
system( paste("time /path/to/software_scripts/bedtools intersect -a ", OUTFILE_5," -b", reference_genome," -wb | grep -v -w chromosome > ", OUTFILE_7_thres_1_genome_v_outliers) );


cat("\nFinished mapping outliers BED files of top 0.05% and 1.0% against reference genome...\n\n");
cat("\nJobs DONE. Bye!\n\n");





