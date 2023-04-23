# After FASTQ analysis, this.
#
# conmmand_line.trimmomatic:
#    java -jar $trimmomatic PE -threads 30 -phred33 -trimlog outfile-trimclip_logfile ${sample}_R1 ${sample}_R2  $paired_dir_path/${sample}_R1_paired  $unpaired_dir_path/${sample}_R1_unpaired  $paired_dir_path/${sample}_R2_paired   $unpaired_dir_path/${sample}_R2_unpaired  ILLUMINACLIP:/node007/share_tools/tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#
# Then, the pipeline below.
# 
# As good practice, check for errors in each BAM file:
# Example:
#    java -Xmx60G -jar $PICARD ValidateSamFile I=$outfiles/${BAM}.paired.dedup_reads.bam  O=$outfiles/${BAM}.paired.dedup_reads.ValidateSamFile.txt   MODE=SUMMARY
#



/node007/share_tools/tools/bwa/bwa mem -M -t 16 \
/node007/users/nlozada/AaegL5/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fasta \
${sample}_1.fastq.gz \
${sample}_2.fastq.gz \
-R "@RG\tID:${sample}\tLB:WholeGenome\tSM:${sample}\tPL:ILLUMINA\tPU:FlowCellId" | \
samtools view -@ 16 -b - > /node007/users/nlozada/AaegL5/alignments/${sample}.bam

/node007/share_tools/tools/java/jre1.7.0_51/bin/java -Xmx32768M -jar /node007/share_tools/tools/picard/AddOrReplaceReadGroups.jar \
INPUT=/node007/users/nlozada/AaegL5/alignments/${sample}.bam \
OUTPUT=/node007/users/nlozada/AaegL5/alignments/${sample}_RRG.bam \
RGID=${sample} \
RGLB=WholeGenome \
RGSM=${sample} \
RGPL=ILLUMINA \
RGPU=FlowCellId

/node007/share_tools/tools/java/jre1.7.0_51/bin/java -Xmx32768M -jar /node007/share_tools/tools/picard/CleanSam.jar \
INPUT=/node007/users/nlozada/AaegL5/alignments/${sample}_RRG.bam \
OUTPUT=/node007/users/nlozada/AaegL5/alignments/${sample}_CLEANED.bam \
VALIDATION_STRINGENCY=SILENT

/node007/share_tools/tools/java/jre1.7.0_51/bin/java -Xmx32768M -jar /node007/share_tools/tools/picard/SortSam.jar \
INPUT=/node007/users/nlozada/AaegL5/alignments/${sample}_CLEANED.bam \
OUTPUT=/node007/users/nlozada/AaegL5/alignments/${sample}_SORTED.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=True

/node007/share_tools/tools/java/jre1.7.0_51/bin/java -Xmx32768M -jar /node007/share_tools/tools/picard/MarkDuplicates.jar \
INPUT=/node007/users/nlozada/AaegL5/alignments/${sample}_SORTED.bam \
OUTPUT=/node007/users/nlozada/AaegL5/alignments/Mark_Duplicates/${sample}.bam \
METRICS_FILE=/node007/users/nlozada/AaegL5/alignments/Mark_Duplicates/${sample}_MD_metrics.txt \
ASSUME_SORTED=True \
CREATE_INDEX=True
