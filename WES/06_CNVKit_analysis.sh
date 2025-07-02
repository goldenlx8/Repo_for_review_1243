#!/bin/bash

#SBATCH --account=
#SBATCH --job-name='CNVKit_pooled_analysis_revised_run'
#SBATCH --output=logs/CNVKit_pooled_analysis_revised_run.log

#SBATCH --partition=standard
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

source activate ~/.conda/envs/cnvkit

ref=data/references/mm39.fa
access=outputs/CNVKit/access.mm39.bed
bam_dir=outputs/GATK_preprocessing
outputDir=outputs/CNVKit
mkdir -p $outputDir
bait=data/references/S0276129/S0276129_Regions_mm39.bed
samples_manifest=outputs/CNVKit_manifest.txt
refFlat=data/references/refFlat.txt

normal_fat_samples=$(cat $samples_manifest | awk '$2 == "Normal_fat" {print "'$bam_dir'/"$1".bam"}' | uniq)
tumor_samples=$(cat $samples_manifest | awk '$2 == "DDLPS" || $2 == "WDLPS" {print "'$bam_dir'/"$1".bam"}' | uniq)

cnvkit.py access $ref -o $access

cnvkit.py batch $tumor_samples --normal $normal_fat_samples \
    --targets $bait --annotate $refFlat \
    --fasta $ref --access $access \
    --output-reference $outputDir/mm39.cnn --output-dir $outputDir \
    --diagram --scatter