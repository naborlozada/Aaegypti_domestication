## Quality control of samples


### FASTQ analysis

```bash
$FILE = *.fastq.gz
for $FILE in *.fastq.gz; do fastqc $FILE -o  $WORKING_DIRECTORY/ ; done
```

### Qualimap analysis

```bash
# INFILE = *.bam
for INFILE in *.bam; do qualimap bamqc  -bam $INFILE --java-mem-size=5G -c -hm 3 -nw 400 -nt 32 -outfile ${newfile}.stats.qualimap -outformat PDF:HTML; done
```

### Mosdepth

```bash
INFILE = *.bam
# Set local directory to call binaries:
export PATH=$PATH:/home/nabor/programs/mosdepth/

# for each chromosomeCHRM_1 jobs:
for INFILE in *.bam; do mosdepth --by AaegL5.regions.bed --thresholds 1,10,20,30  $WORKING_DIRECTORY/${newfile}.depth $INFILE; done
```

