#!/bin/bash

#SBATCH --account=
#SBATCH --job-name='merging_bams'
#SBATCH --output=logs/bwa_%a.log

#SBATCH --partition=standard
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

#SBATCH --array=1-23

## Loading Important Modules from U-M Greatlakes HPC
ml samtools

## Setting up variables for files and directories
samples_manifest=outputs/bwa_manifest.txt
sample=$(cat $samples_manifest | cut -f1 | uniq | sed -n $[SLURM_ARRAY_TASK_ID]p)
outputDir=outputs/GATK_preprocessing
mkdir -p ${outputDir}

## Merge bams of different samples
samtools merge -r -o ${outputDir}/${sample}.bam ${outputDir}/${sample}*

## Removing old files
rm -vf ${outputDir}/${sample}*_BQSR.bam

samtools index ${outputDir}/${sample}.bam