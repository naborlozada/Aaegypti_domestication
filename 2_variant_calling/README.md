## Variant Calling

This directory contains one single bash script file:
* **`variant_calling_and_hardfiltering.sh`:** Bash script used to perform 1) a variant calling in a single population, producing a single VCF file with all raw-unfiltered variants (indels/SNPs), 2) removes indels with a close proximity to SNPs, 3) extract SNPs from raw filtered variants, 4) filters SNPs, and 5) extract biallelic SNPs.  

This bash script is adapted to be used in the slurm jobs cluster system to run several jobs with multiple tasks in parallel. All tasks (1-5) run one-by-one in that specific order.
