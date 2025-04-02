#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RAISD example command:
#    RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -n uganda_kichwamba.nr.chrm1.biallelic.geno80.win50.output -I uganda.kichwamba.nr.chrm1.biallelic.geno80.vcf.gz
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


rm(list=ls());

.libPaths("/scr/node007/nlozada/software_scripts/R/v3.6.2/lib64/R/library/");
#.libPaths();



Sys.time();

suppressPackageStartupMessages(require(qvalue));
suppressPackageStartupMessages(library(qqman));
suppressPackageStartupMessages(library(dplyr));
#suppressPackageStartupMessages(library(tibble));




# PARSE INFILE
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Report   <- "RAiSD_Report.uganda_kichwamba.nr.chrm1.biallelic.geno80.win50.output.1";

# chromosome: 
chromTAG <- as.numeric(1);


# parse info
d       <- read.table(paste(Report[1]), header=F, skip=0)[,1:7];
data    <- data.frame(pos=c(), value=c(), chr=c());
tmp.dat <- data.frame(pos=d[,1], value=d[,7]*1, chr=rep(chromTAG, length(d[,1])));
data    <- rbind(data, tmp.dat);

#head(data);

# parse data
snp <- 1:dim(data)[1];
mydf <- data.frame(snp, data);
head(mydf);





# OUTLIERS DETECTION   
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


raisd_full_site_report.outliers.mahalanobis = data.frame(mu_stat = mydf$value);

##  *** Mahalanobis test ***
raisd_full_site_report.outliers.mahalanobis$mahalanobis <- stats::mahalanobis(raisd_full_site_report.outliers.mahalanobis, base::colMeans(raisd_full_site_report.outliers.mahalanobis), stats::cov(raisd_full_site_report.outliers.mahalanobis));
# get pvalues
raisd_full_site_report.outliers.mahalanobis$pvalue  <- pchisq(raisd_full_site_report.outliers.mahalanobis$mahalanobis, df=1, lower.tail=FALSE);

# Merge columns of both dataframes
raisd_full_site_report.mahalanobis.full_info <- cbind(mydf, raisd_full_site_report.outliers.mahalanobis);




## *** QVALUES ***

# transform pvalues:
qvalues <- qvalue(raisd_full_site_report.mahalanobis.full_info$pvalue)$qvalues;

qvalues_df <- as.data.frame(qvalues);

# alpha threshold:
alpha <- 0.05;

outliers.qvalues <- which(qvalues < alpha);
length(outliers.qvalues);

# merge all info:
myData_stats <- cbind(raisd_full_site_report.mahalanobis.full_info,qvalues_df);

# define outliers
myData_stats_qvals_categories <- myData_stats %>% mutate(outlier_qvalFDR5 = ifelse(qvalues < alpha, "outlier", "background"));


write.table(myData_stats_qvals_categories, file=myData_stats_all_signif_categories.TXT, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE);

