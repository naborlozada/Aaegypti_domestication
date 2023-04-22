# Realignment and recalibration of BAM files:


Infiles and directories

```bash
export TMPDIR=/node007/users/nabor/TMP_DIR/

tmpdir=/node007/users/nabor/TMP_DIR/
GATK=/node007/users/nabor/programs/gatk
reference_genome=/node007/users/nabor/nabor/aedes_aegypti/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta
WORKING_DIR=/node007/users/nabor/nabor/aedes_aegypti/
varDB_indels=/node007/users/nabor/nabor/aedes_aegypti/indels.Db.vcf.gz
varDB_snps=/node007/users/nabor/nabor/aedes_aegypti/snps.Db.vcf.gz

```


## Indels

Step 1:

```bash
# $FILE= your deduplicated BAM file  // OR in a cluster:  FILE=${FILES[$SLURM_ARRAY_TASK_ID]} 
# ${newfile}= new name for your BAM file

time java -XX:ParallelGCThreads=32 -Xmx40G -jar $GATK \
    -T RealignerTargetCreator \
    -R $reference_genome \
    -I $FILE \
    -known $knownINDELs \
    -o $WORKING_DIR/${newfile}.paired.dedup_reads.indels.intervals \
    --num_threads 32
```

Step 2:

```bash
time java -Djava.io.tmpdir=${tmpdir} -XX:ParallelGCThreads=1 -Xmx5G -jar $GATK \
    -T IndelRealigner \
    -R $reference_genome \
    -I $FILE \
    -known $knownINDELs \
    -targetIntervals $infiles/${newfile}.paired.dedup_reads.indels.intervals \
    -o $WORKING_DIR/${newfile}.paired.dedup_reads.indels.realigned.bam

wait;

# time $samtools index -@ 16 $WORKING_DIR/${newfile}.paired.dedup_reads.indels.realigned.bam
```

Step 3:

```bash
time java -Djava.io.tmpdir=${tmpdir} -Xmx40G -XX:+UseParallelGC -XX:ParallelGCThreads=16 -jar $picard FixMateInformation \
    INPUT=$FILE \
    OUTPUT=$WORKING_DIR/${newfile}.paired.dedup_reads.indels.realigned.fixed_mate.bam \
    ADD_MATE_CIGAR=true \
    TMP_DIR=$tmpdir

wait;

# time $samtools index -@ 16 $WORKING_DIR/${newfile}.paired.dedup_reads.indels.realigned.bam
```


## Recalibration

```bash
# FIRST PASS:
time $GATK \
    --java-options "-Xmx40G -Djava.io.tmpdir=${tmpdir} -XX:ParallelGCThreads=3 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" BaseRecalibrator \
    --input $FILE \
    --reference $reference_genome \
    --known-sites $varDB_indels \
    --known-sites $varDB_snps \
    --output $WORKING_DIR/${newfile}.dedup.indels.realn.fixed_mate.bqsr1.before.tbl \
    --tmp-dir $tmpdir

wait;


time $GATK \
    --java-options "-Xmx40G -Djava.io.tmpdir=${tmpdir} -XX:ParallelGCThreads=3 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyBQSR \
    --reference $reference_genome \
    --input $FILE \
    --bqsr-recal-file $WORKING_DIR/${newfile}.dedup.indels.realn.fixed_mate.bqsr1.before.tbl \
    --output $WORKING_DIR/${newfile}.dedup.indels_realn.bqsr.recal.bam \
    --tmp-dir $tmpdir

wait;
```

```bash
# SECOND PASS:
time $GATK \
    --java-options "-Xmx40G -Djava.io.tmpdir=${tmpdir} -XX:ParallelGCThreads=3 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" BaseRecalibrator \
    --input $WORKING_DIR/${newfile}.dedup.indels_realn.bqsr.recal.bam \
    --reference $reference_genome \
    --known-sites $varDB_indels \
    --known-sites $varDB_snps \
    --output $WORKING_DIR/${newfile}.dedup.indels.realn.fixed_mate.bqsr2.after.tbl \
    --tmp-dir $tmpdir

wait;

# If you reach this point (after a lot of effort), then you got your alignments recalibrated!


# Compare the First vs Second round:
time $GATK \
    --java-options "-Xmx40G -Djava.io.tmpdir=${tmpdir} -XX:ParallelGCThreads=3 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" AnalyzeCovariates \
    -before $WORKING_DIR/${newfile}.dedup.indels.realn.fixed_mate.bqsr1.before.tbl \
    -after $WORKING_DIR/${newfile}.dedup.indels.realn.fixed_mate.bqsr2.after.tbl \
    -plots $WORKING_DIR/${newfile}.recalibration_plots.bqsr.before_after.pdf \
    -csv $WORKING_DIR/${newfile}.recalibration_plots.bqsr.before_after.csv

wait;

# ALmost done... Maybe one last thing: 
# Check any potential error in your alignment

time java -Djava.io.tmpdir=${tmpdir} -Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=16 -jar $picard ValidateSamFile \
    INPUT=$FILE \
    OUTPUT=$WORKING_DIR/${newfile}.dedup.indels_realn.bqsr.recal.ValidateSamFile.txt \
    MODE=SUMMARY \
    TMP_DIR=$tmpdir

wait;

```


