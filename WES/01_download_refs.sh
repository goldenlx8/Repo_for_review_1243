########################################
### Download References Genome from UCSC
########################################
mkdir -p data/references
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz -O data/references/mm39.fa.gz
gunzip data/references/mm39.fa.gz

###########################
### Index References Genome
###########################
ref=data/references/mm39.fa
bwa index $ref
samtools faidx $ref
gatk CreateSequenceDictionary -R $ref
## to run the last line you want to install cnvkit
## and activate its conda env

###########################
### refFlat mm39
###########################
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz -o data/references/refFlat.txt.gz
gunzip data/references/refFlat.txt.gz

################################
## Download Known sites for BQSR
################################
## mm10

# ## for mm10 as described in (https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md)
# ## NCBI SNPs database
# vcf=data/references/00-All.vcf.gz
# wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz -O $vcf
# vcf_new=data/references/00-All.vcf_mm10.gz
# zcat $vcf | sed 's/^\([0-9XY]\)/chr\1/' > $vcf_new
# rm -fv $vcf
# ## MGP database for indels
# wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz \
# -O data/references/mgp.v5.indels.vcf.gz
# # Filter for passing variants with chr added:
# # adjust header
# zcat data/references/mgp.v5.indels.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 \
# | grep -v "#contig" | grep -v "#source" \
# > data/references/mgp.v5.indels.pass.chr.vcf
# # keep only passing and adjust chromosome name
# zcat data/references/mgp.v5.indels.vcf.gz | grep -v "^#" | cut -f 1-8 \
# | grep -w "PASS" | sed 's/^\([0-9MXY]\)/chr\1/' \
# >> data/references/mgp.v5.indels.pass.chr.vcf
# # Sort VCF (automatically generated index has to be deleted due to a known bug):
# java -Xms16G -Xmx16G -jar $PICARDLIB/picard.jar SortVcf VERBOSITY=WARNING \
# SD=genome.dict \
# I=data/references/mgp.v5.indels.pass.chr.vcf \
# O=data/references/mgp.v5.indels.pass.chr.sort.vcf
# rm -fv data/references/mgp.v5.indels.pass.chr.sort.vcf.idx

## mm39

## Since we are using mm39 here, and the NCBI database is based on mm10,
## and the MGP website is archived, we will use the datasets from the
## https://www.mousegenomes.org/ for mm39

## Mouse Genome Project SNPs database
mpg_snps=data/references/mgp_REL2021_snps.vcf.gz
mod_mpg_snps=data/references/mgp_REL2021_snps.chr.vcf
final_mpg_snps=data/references/mgp_REL2021_snps.chr.sorted.vcf
wget ftp://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz \
-O $mpg_snps
# adjust header
zcat $mpg_snps | head -1000 | grep "^#" | cut -f 1-8 \
| grep -v "#contig" | grep -v "#source" \
> $mod_mpg_snps
# keep only passing and adjust chromosome name
zcat $mpg_snps | grep -v "^#" | cut -f 1-8 \
| grep -w "PASS" | sed 's/^\([0-9MXY]\)/chr\1/' \
>> $mod_mpg_snps
# Sort VCF
java -Xms16G -Xmx16G -jar $PICARDLIB/picard.jar SortVcf \
SD=data/references/mm39.dict \
I=$mod_mpg_snps \
O=$final_mpg_snps

rm -vf $mpg_snps $mod_mpg_snps

# Mouse Genome Project INDEL database
mpg_indels=data/references/mgp_REL2021_indels.vcf.gz
mod_mpg_indels=data/references/mgp_REL2021_indels.chr.vcf
final_mpg_indels=data/references/mgp_REL2021_indels.chr.sorted.vcf
wget ftp://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_indels.vcf.gz \
-O $mpg_indels
# adjust header
zcat $mpg_indels | head -1000 | grep "^#" | cut -f 1-8 \
| grep -v "#contig" | grep -v "#source" \
> $mod_mpg_indels
# keep only passing and adjust chromosome name
zcat $mpg_indels | grep -v "^#" | cut -f 1-8 \
| sed 's/^\([0-9MXY]\)/chr\1/' \
>> $mod_mpg_indels
# Sort VCF
java -Xms16G -Xmx16G -jar $PICARDLIB/picard.jar SortVcf \
SD=data/references/mm39.dict \
I=$mod_mpg_indels \
O=$final_mpg_indels

rm -vf $mpg_indels $mod_mpg_indels

###########################
### WES Capture bed file
###########################
## WES capture kit bed file is downloaded as follows
## (https://kb.10xgenomics.com/hc/en-us/articles/115004150923-Where-can-I-find-the-Agilent-Target-BED-files-)
## Since the bait file is based on mm9, liftover tool (https://genome.ucsc.edu/cgi-bin/hgLiftOver) was used to
## map the coordinates to mm39

wget --timestamping \
    'ftp://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm39.over.chain.gz' \
    -O data/references/mm9ToMm39.over.chain.gz

gunzip data/references/mm9ToMm39.over.chain.gz

tail -n +3 data/references/S0276129/S0276129_Regions.bed > data/references/S0276129/S0276129_Regions_no_header.bed

./liftOver \
data/references/S0276129/S0276129_Regions_no_header.bed \
data/references/mm9ToMm39.over.chain \
data/references/S0276129/S0276129_Regions_mm39.bed unMapped