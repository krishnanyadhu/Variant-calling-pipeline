#!/bin/bash

set -o errexit
set -o nounset

for i in $(cat data.txt)
do

#make main output direcotry
mkdir -p output/${i}

mkdir -p output/${i}/reads
#download_data_fastq_format in to the output/reads directory
#if it is a binary file use fasterq-dump to covert them into reads
fasterq-dump data/${i} -O output/${i}/reads -e 30 --split-files

#quality_check_fastqc
#make_output_directory for fastqc output
mkdir -p output/${i}/fastqc_output
fastqc -t 30 output/${i}/reads/${i}_1.fastq -o output/${i}/fastqc_output
fastqc -t 30 output/${i}/reads/${i}_2.fastq -o output/${i}/fastqc_output

#quality control and adapter trimmming using fastp
mkdir -p output/${i}/trimmed_reads
fastp -i output/${i}/reads/${i}_1.fastq -I output/${i}/reads/${i}_2.fastq -o output/${i}/trimmed_reads/${i}_1_trimmed.fq -O output/${i}/trimmed_reads/${i}_2_trimmed.fq -j output/${i}/trimmed_reads/${i}_fastp.json --detect_adapter_for_pe -h output/${i}/trimmed_reads/${i}_fastp.html -w 16

#quality check after quality and adapter trimming
mkdir -p output/${i}/fastp_fastqc_output
fastqc -t 30 output/${i}/trimmed_reads/${i}_1_trimmed.fq -o output/${i}/fastp_fastqc_output
fastqc -t 30 output/${i}/trimmed_reads/${i}_2_trimmed.fq -o output/${i}/fastp_fastqc_output


#align with indexed reference genome
mkdir -p  output/${i}/BWA_alignment
bwa mem -t 30 Ref/hg38.fa output/${i}/trimmed_reads/${i}_1_trimmed.fq output/${i}/trimmed_reads/${i}_2_trimmed.fq > output/${i}/BWA_alignment/${i}_alignment.sam

#SAMtoBAM
mkdir -p output/${i}/BWA_alignment_conversion
samtools view -@30 -bS output/${i}/BWA_alignment/${i}_alignment.sam -o output/${i}/BWA_alignment_conversion/${i}_alignment.bam

#sortbam
mkdir -p output/${i}/BWA_alignment_sort
samtools sort -@30 output/${i}/BWA_alignment_conversion/${i}_alignment.bam -o  output/${i}/BWA_alignment_sort/${i}_alignment_sort.bam

#add_read_group
mkdir -p output/${i}/add_read_group
java -Xmx100g -jar tools/picard-tools/picard.jar AddOrReplaceReadGroups \
INPUT=output/${i}/BWA_alignment_sort/${i}_alignment_sort.bam \
OUTPUT=output/${i}/add_read_group/${i}_alignment_RG.bam \
SORT_ORDER=coordinate \
RGID=${i} \
RGLB=${i} \
RGPL=illumina \
RGPU=${i} \
RGSM=${i} \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT

#remove_PCR_duplicates
mkdir -p output/${i}/PCR_duplicates_removed
java -Xmx100g -jar tools/picard-tools/picard.jar MarkDuplicates \
INPUT=output/${i}/add_read_group/${i}_alignment_RG.bam \
OUTPUT=output/${i}/PCR_duplicates_removed/${i}_alignment_PCR.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE=output/${i}/PCR_duplicates_removed/${i}_alignment_PR.Metrics \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true

#BQSR
#baserecalibrator
mkdir -p output/${i}/BQSR
gatk BaseRecalibrator -I output/${i}/PCR_duplicates_removed/${i}_alignment_PCR.bam \
-R Ref/hg38.fa \
--known-sites Ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--known-sites Ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O output/${i}/BQSR/${i}_recal_data.table

#ApplyBQSR
gatk ApplyBQSR \
-R Ref/hg38.fa \
-I output/${i}/PCR_duplicates_removed/${i}_alignment_PCR.bam \
--bqsr-recal-file output/${i}/BQSR/${i}_recal_data.table \
-O output/${i}/BQSR/${i}_recaliberated.bam


#Extract chromosome names
chromosomes=$(samtools idxstats output/${i}/PCR_duplicates_removed/${i}_alignment_PCR.bam | cut -f1 | grep -E '^chr[0-9XYM]+$')
#chromosomes=$(samtools idxstats output/${i}/PCR_duplicates_removed/${i}_alignment_PCR.bam | cut -f1 | sort | uniq )
echo $chromosomes

#parallel_variant_calling
mkdir -p output/${i}/TMP_DIR
echo "$chromosomes" | parallel -j 10 "gatk PrintReads --input output/${i}/BQSR/${i}_recaliberated.bam -L "{}" --output output/"${i}"/TMP_DIR/{}_only.bam && gatk HaplotypeCaller -R Ref/hg38.fa -I output/"${i}"/TMP_DIR/"{}"_only.bam -L "{}" -O output/"${i}"/TMP_DIR/"{}"_GATK.vcf.gz --native-pair-hmm-threads 4"


#Merge all chromosome VCFs
mkdir -p output/${i}/variant_calling
VCF_FILES=$(ls output/${i}/TMP_DIR/chr*_GATK.vcf.gz | awk '{print "-I " $1}')
java -Xmx100g -jar tools/picard-tools/picard.jar MergeVcfs -O output/${i}/variant_calling/${i}_GATK.vcf.gz $VCF_FILES

rm -r output/${i}/TMP_DIR
done
