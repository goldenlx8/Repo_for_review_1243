#!/bin/bash

#SBATCH --account=
#SBATCH --job-name='bwa'
#SBATCH --output=logs/bwa_%a.log

#SBATCH --partition=standard
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

#SBATCH --array=1-47

## Loading Important Modules from U-M Greatlakes HPC
ml samtools
ml bwa
ml gatk
ml picard-tools

## Setting up variables for files and directories
samples_manifest=outputs/GATK_preprocess_manifest.txt
ref=data/references/mm39.fa
outputDir=outputs/GATK_preprocessing
mkdir -p ${outputDir}

## Subsetting a sample
sample=$(cat $samples_manifest | cut -f2 | sed -n $[SLURM_ARRAY_TASK_ID]p)
R1=$(cat $samples_manifest | cut -f3 | sed -n $[SLURM_ARRAY_TASK_ID]p)
R2=$(cat $samples_manifest | cut -f4 | sed -n $[SLURM_ARRAY_TASK_ID]p)
known_sites1=data/references/mgp_REL2021_indels.chr.sorted.vcf
known_sites2=data/references/mgp_REL2021_snps.chr.sorted.vcf

## Logging information
echo "::> Analyzing Sample $sample"
echo -e "::> R1 path $R1"
echo -e "::> R2 path $R2"
echo -e "::> Results will be stored in ${outputDir}"

## Alignment using BWA
echo "::> Aligning sample using BWA <::"
bwa mem $ref $R1 $R2 -t 16  -o ${outputDir}/${sample}.sam
samtools view -b --threads 16 --verbosity 1 ${outputDir}/${sample}.sam -o ${outputDir}/${sample}.bam
rm -vf ${outputDir}/${sample}.sam

## GATK Preprocessing
echo "::> GATK Preprocessing <::"
gatk AddOrReplaceReadGroups -I ${outputDir}/${sample}.bam \
-O ${outputDir}/${sample}_RG.bam \
-LB ${sample} \
-PL Illumina \
-PU ${sample} \
-SM ${sample} \
--VERBOSITY ERROR
rm -vf ${outputDir}/${sample}.bam

samtools sort -O bam -o ${outputDir}/${sample}_sorted.bam ${outputDir}/${sample}_RG.bam
rm -vf ${outputDir}/${sample}_RG.bam

samtools index -b ${outputDir}/${sample}_sorted.bam

gatk MarkDuplicates -I ${outputDir}/${sample}_sorted.bam \
-O ${outputDir}/${sample}_MD.bam \
-M ${outputDir}/${sample}_MD.txt \
--VERBOSITY ERROR
rm -vf ${outputDir}/${sample}_sorted.bam*

gatk BaseRecalibrator \
-I ${outputDir}/${sample}_MD.bam \
-R $ref \
--known-sites $known_sites1 \
--known-sites $known_sites2 \
-O ${outputDir}/${sample}_DupStat.table \
--verbosity ERROR

gatk ApplyBQSR -R $ref \
-I ${outputDir}/${sample}_MD.bam \
--bqsr-recal-file ${outputDir}/${sample}_DupStat.table \
-O ${outputDir}/${sample}_BQSR.bam \
--verbosity ERROR

rm -vf ${outputDir}/${sample}_MD.*
rm -vf ${outputDir}/${sample}_DupStat.table