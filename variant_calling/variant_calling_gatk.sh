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
#SBATCH --time=
#SBATCH --error=preprocessing_reads.%A_%a.stderr.log
#SBATCH --output=preprocessing_reads.%A_%a.stderr.log




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



starttime=`date +"%s"`

export TMPDIR=/node007/users/nabor/TMP_DIR/

tmpdir=/node007/users/nabor/TMP_DIR/
GATK=/node007/users/nabor/programs/gatk
reference_genome=/node007/users/nabor/nabor/aedes_aegypti/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta
outfiles=/node007/users/nabor/nabor/aedes_aegypti/
varDB_indels=/node007/users/nabor/nabor/aedes_aegypti/indels.Db.vcf.gz
varDB_snps=/node007/users/nabor/nabor/aedes_aegypti/snps.Db.vcf.gz



echo
echo JOB_START: `date`
echo
echo " *** Variant calling using GATK (variant calling and filering): ***"
echo
echo "##############################################################################################################################"
echo
echo

# Include all BAM samples or a single BAM file
time java -Djava.io.tmpdir=${tmpdir} -Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -jar $GATK \
    -T UnifiedGenotyper \
    --reference_sequence $refGENOME \
    --input_file $WORKING_DIRECTORY/BAM_sample1.recal.bam \
    --input_file $WORKING_DIRECTORY/BAM_sample2.recal.bam \
    --input_file $WORKING_DIRECTORY/BAM_sample3.recal.bam \
    --dbsnp $varDB_snps \
    --genotype_likelihoods_model BOTH \
    --downsample_to_coverage 200 \
    --min_base_quality_score 20 \
    --num_threads 1 \
    -o $WORKING_DIRECTORY/mySNPs.raw.vcf

wait;

# -----------------------------------------------------------------------------
# SNP extraction.
time java -Djava.io.tmpdir=${tmpdir} -Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -jar $GATK \
    -T SelectVariants \
    -R $refGENOME \
    -V $WORKING_DIRECTORY/mySNPs.raw.vcf.gz \
    -selectType SNP \
    -o $WORKING_DIRECTORY/mySNPs.raw.snps.vcf
wait;


bgzip -c $WORKING_DIRECTORY/Population_Name.raw.snps.vcf > $WORKING_DIRECTORY/mySNPs.raw.snps.vcf.gz
wait;

bcftools index -t $WORKING_DIRECTORY/mySNPs.raw.snps.vcf.gz
wait;

# -----------------------------------------------------------------------------
# Remove SNPs with a close proximity to indes of 10 bp:
bcftools filter --threads 3 --SnpGap 10 -Oz -o $WORKING_DIRECTORY/mySNPs.output1.vcf.gz  $WORKING_DIRECTORY/mySNPs.raw.snps.vcf.gz
wait;


# Hard filtering options for SNPs using GATK:
time java -Djava.io.tmpdir=${tmpdir} -Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -jar $GATK \
    -T VariantFiltration \
    -R $refGENOME \
    -V $WORKING_DIRECTORY/mySNPs.output1.vcf.gz \
    --filterExpression "QUAL < 19.9 || QD < 2.0 || FS > 60.0 || MQ < 19.9 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \ 
    --filterName "my_snp_filter" \ 
    -o $WORKING_DIRECTORY/mySNPs.output2.vcf.gz
wait;

bcftools view -f PASS  -Oz -o $WORKING_DIRECTORY/mySNPs.output3.vcf.gz  $WORKING_DIRECTORY/mySNPs.output2.vcf.gz
wait;

# -----------------------------------------------------------------------------
# Extraction of biallelic snps
bcftools view -m2 -M2 -Oz -o mySNPs.biallelic.vcf.gz  $WORKING_DIRECTORY/mySNPs.output3.vcf.gz
wait;

# remove tmp files
rm mySNPs.output*.vcf*
wait;


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
