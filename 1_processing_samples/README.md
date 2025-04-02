## Processing samples

This directory contains three main files:
* **`quality_control.md`:** Commands used to performed a quality control in reads (\*.fastq.gz) and alignments (BAM format) across the whole genome and to specific genomic regions (cds, introns, transcripts, etc.). These commands run one single task at time in a loop procees using the function `for`.


* **`alignments.sh`:** Bash script used to 1) trim adapters, 2) align paired ended reads to the reference genome of AaegL5, 3) sort aligned reads, and 4) mark duplicated reads (deduping process). This script is specifically adapted to be used it in a slurm jobs cluster system, and all these tasks (1-4) run one-by-one in that specific numeric order.


* **`recalibration.md`:** Set of GATK command line functions used to 1) re-align indels and 2) recalibrate alignments in BAM format. These commands 
  
