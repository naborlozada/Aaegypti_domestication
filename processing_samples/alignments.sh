#!/usr/bin/bash

# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada Chavez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #


#SBATCH --account=
#SBATCH --partition=
#SBATCH --cpus-per-task=32
#SBATCH --mem=
#SBATCH --ntasks=1
#SBATCH --array=0-10

#SBATCH --job-name=
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

#SBATCH --time=7-00:00:00
#SBATCH --error=preprocessing_reads.%A_%a.stderr.log
#SBATCH --output=preprocessing_reads.%A_%a.stderr.log


# #SBATCH --ntasks-per-node=64

echo
echo
echo "            **** SBATCH JOBS SUBMISSION wih SLURM system ****                 "
echo 
echo

echo "------------------------------------------------------------------------------"
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "JOB submission name            = $SLURM_JOB_NAME"
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "User account name              = $SLURM_JOB_ACCOUNT"
echo "Partition name                 = $SLURM_JOB_PARTITION"
echo "Export environment DIR         = $SLURM_EXPORT_ENV"
echo "Task PID                       = $SLURM_TASK_PID $SLURM_JOB_ID"
echo "------------------------------------------------------------------------------"
echo
echo
echo




export TMPDIR=node007/users/nabor/TMP_DIR/

# // paths & programs //
# --------------------------------------------------
# programs
trimmomatic=/node007/users/nabor/tools/programs/bioinformatics/Trimmomatic-0.39/trimmomatic-0.39.jar
SAMTOOLS=/node007/users/nabor/tools/programs/bioinformatics/samtools-1.10/bin/samtools
PICARD=/node007/users/nabor/tools/programs/bioinformatics/picard.jar
BWA=/node007/users/nabor/tools/programs/bioinformatics/bwa-0.7.17/bwa

# infiles
reference_genome=Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta
reference_dir=/node007/users/nabor/tools/aedes_aegypti/genome/reference_genome
reads_path=/node007/users/nabor/tools/aedes_aegypti/results/alignments/POPULATIONS/reads.raw_n_sort
tmp_files=/node007/share_tools/tools/TMP_DIR/

# outfiles
paired_dir_path=/node007/users/nabor/tools/aedes_aegypti/results/alignments/POPULATIONS/reads.raw_n_sort
unpaired_dir_path=/node007/users/nabor/tools/aedes_aegypti/results/alignments/POPULATIONS/reads.raw_n_sort
sortedDIR=/node007/users/nabor/tools/aedes_aegypti/results/alignments/POPULATIONS/reads.raw_n_sort
dedupDIR=/node007/users/nabor/tools/aedes_aegypti/results/alignments/POPULATIONS/bam.dedups

# infiles
FILES1=(/node007/users/nabor/tools/aedes_aegypti/results/alignments/POPULATIONS/reads.raw_n_sort/*_1.fastq.gz)
FILE1=${FILES1[$SLURM_ARRAY_TASK_ID]}
FILE2=$(basename ${FILE1} _1.fastq.gz)_2.fastq.gz

newfile="$(basename $FILE1 _1.fastq.gz)"


starttime=`date +"%s"`


echo
echo JOB_START: `date`
echo
echo " *** Processing samples (trimming,align,sort,dedup,validate): $FILE1 $FILE2 ***"
echo
echo "##############################################################################################################################"
echo
echo
echo "READS FILES R1 R2: $FILE1 $FILE2 ..."
echo 
# -----------------------------------------------------------------------------
echo "[1] // Trimmomatic: Trimming reads adapters //"
echo "command_line: java -Djava.io.tmpdir=$tmp_files -jar $trimmomatic PE -threads 32 -phred33  -trimlog $sortedDIR/${newfile}.R1_R2.trimmomatic.logfile.txt  $FILE1  $reads_path/$FILE2  $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz  $unpaired_dir_path/${newfile}.R1.unpaired.trimclip.fastq.gz  $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz $unpaired_dir_path/${newfile}.R2.unpaired.trimclip.fastq.gz ILLUMINACLIP:/home/chavez/nabor/programs/bioinformatics/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
echo

time java -Djava.io.tmpdir=$tmp_files -jar $trimmomatic PE -threads 32 -phred33 \
            -trimlog $sortedDIR/${newfile}.R1_R2.trimmomatic.logfile.txt \
            $FILE1 \
            $reads_path/$FILE2 \
            $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz \
            $unpaired_dir_path/${newfile}.R1.unpaired.trimclip.fastq.gz \
            $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz \
            $unpaired_dir_path/${newfile}.R2.unpaired.trimclip.fastq.gz \
            ILLUMINACLIP:/home/chavez/nabor/programs/bioinformatics/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
wait;
# -----------------------------------------------------------------------------
echo
echo "[2] // BWA: Alignment of reads to  he reference genome // "
echo "command_line: $BWA mem -t 32  $reference_dir/$reference_genome  $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz  $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz  -R '@RG\tID:'$newfile'\tSM:'$newfile'\tPL:illumina' | $SAMTOOLS view -b --threads 32 - > $sortedDIR/${newfile}.paired.aligned_raw_reads.bam "
echo

time $BWA mem -t 32 \
            $reference_dir/$reference_genome \
            $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz \
            $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz \
            -R '@RG\tID:'$newfile'\tSM:'$newfile'\tPL:illumina' | $SAMTOOLS view -b --threads 32 - > $sortedDIR/${newfile}.paired.aligned_raw_reads.bam
wait;
# -----------------------------------------------------------------------------
echo
echo "[3] // SortSam: Sort alignment // "
echo "command_line: java -Djava.io.tmpdir=$tmp_files -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD SortSam INPUT=$sortedDIR/${newfile}.paired.aligned_raw_reads.bam OUTPUT=$sortedDIR/${newfile}.paired.aligned_reads_sorted.bam SORT_ORDER=coordinate  TMP_DIR=$tmp_files "
echo

time java -Djava.io.tmpdir=$tmp_files -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD SortSam \
            INPUT=$sortedDIR/${newfile}.paired.aligned_raw_reads.bam \
            OUTPUT=$sortedDIR/${newfile}.paired.aligned_reads_sorted.bam \
            SORT_ORDER=coordinate \
            TMP_DIR=$tmp_files
wait;
# -----------------------------------------------------------------------------
echo
echo "[4] // MarkDuplicates //"
echo "command_line:  java -Djava.io.tmpdir=$tmp_files -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates INPUT=$sortedDIR/${newfile}.paired.aligned_reads_sorted.bam OUTPUT=$dedupDIR/${newfile}.paired.dedup_reads.bam METRICS_FILE=$dedupDIR/${newfile}.paired.dedup_reads.metrics.txt TMP_DIR=$tmp_files "
echo

time java -Djava.io.tmpdir=$tmp_files -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates \
            INPUT=$sortedDIR/${newfile}.paired.aligned_reads_sorted.bam \
            OUTPUT=$dedupDIR/${newfile}.paired.dedup_reads.bam \
            METRICS_FILE=$dedupDIR/${newfile}.paired.dedup_reads.metrics.txt \
            TMP_DIR=$tmp_files
wait;
# -----------------------------------------------------------------------------
echo
echo "[5] // Make BWA-based index... //"
echo "command_line: $BWA index $dedupDIR/${newfile}.paired.dedup_reads.bam "
echo

$BWA index $dedupDIR/${newfile}.paired.dedup_reads.bam
wait;
# -----------------------------------------------------------------------------
echo
echo "[6] // Make samtools-based index... //"
echo "command_line: $SAMTOOLS index -@ 64 $dedupDIR/${newfile}.paired.dedup_reads.bam"
echo

$SAMTOOLS index -@ 64 $dedupDIR/${newfile}.paired.dedup_reads.bam
wait;
# -----------------------------------------------------------------------------
echo
echo "[7] // Validate SAM/BAM error in alignment... //"
echo "command_line: java -Djava.io.tmpdir=$tmp_files -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD ValidateSamFile I=$dedupDIR/${newfile}.paired.dedup_reads.bam O=$dedupDIR/${newfile}.paired.dedup_reads.ValidateSamFile.txt MODE=SUMMARY TMP_DIR=$tmp_files "
echo

time java -Djava.io.tmpdir=$tmp_files -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD ValidateSamFile \
            I=$dedupDIR/${newfile}.paired.dedup_reads.bam \
            O=$dedupDIR/${newfile}.paired.dedup_reads.ValidateSamFile.txt \
            MODE=SUMMARY \
            TMP_DIR=$tmp_files
wait;
# -----------------------------------------------------------------------------

endtime=`date +"%s"`
duration=$((endtime - starttime))

echo
echo
echo "##############################################################################################################################"
echo
echo JOB_ENDED: `date`
echo
echo
echo "STAT:startTime:$starttime"
echo "STAT:doneTime:$endtime"
echo "STAT:runtime:$duration"
echo
echo



