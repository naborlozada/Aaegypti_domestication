## Quality control of samples


### FASTQ analysis

Check the quality of the reads ion the FASTQ files before the alignment them against the reference genome.

```bash
$FILE = *.fastq.gz   
for $FILE in *.fastq.gz; do fastqc $FILE -o  $WORKING_DIRECTORY/ ; done

## or replace: *.fastq.gz   for this 
##  /full_dir_path/fastqs/*.fastq.gz;
```

### Qualimap analysis

Check the quality of alignment after aligning the reads against the reference genome.

```bash
## INFILE = *.bam
for INFILE in *.bam; do qualimap bamqc  -bam $INFILE --java-mem-size=5G -c -hm 3 -nw 400 -nt 32 -outfile ${newfile}.stats.qualimap -outformat PDF:HTML; done

## or replace: *.bam   for this 
##  /full_dir_path/bams/*.bam;
```

### Mosdepth

Check the coverage of specific genomic regions (genes, transcripts, cds, exons) after the reads alignment against the reference genome.

```bash
INFILE = *.bam
# Set local directory to call binaries:
export PATH=$PATH:/home/nabor/programs/mosdepth/

# for each chromosomeCHRM_1 jobs:
for INFILE in *.bam; do mosdepth --by AaegL5.regions.bed --thresholds 1,10,20,30  $WORKING_DIRECTORY/${newfile}.depth $INFILE; done
```

