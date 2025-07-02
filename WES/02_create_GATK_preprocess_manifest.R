setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(tidyverse)
library(janitor)

dir.create("outputs", showWarnings = F, recursive = T)

fastq_manifest <- read.table("data/fastq_manifest.txt")
fastq_manifest <- lapply(fastq_manifest$V2 , function(x){
  paths <- list.files(x, full.names = T, pattern = '.fq.gz')
  sub_sample <- list.files(x, pattern = '.fq.gz')
  return(data.frame(sub_sample, paths))
}) %>%
  `names<-`(fastq_manifest$V1) %>%
  bind_rows(.id = 'sample_id')
fastq_manifest$read <- str_extract(fastq_manifest$sub_sample, "[1-2].fq.gz$",) %>%
  gsub(".fq.gz","", .)
fastq_manifest$sub_sample <- gsub("_1.fq.gz|_2.fq.gz","", fastq_manifest$sub_sample)
fastq_manifest <- reshape(fastq_manifest, idvar=c("sample_id", "sub_sample"),
                          timevar = "read", direction="wide") %>%
  `colnames<-`(c('sample_id', 'sub_sample', 'R1', 'R2'))

samples_manifest <- readxl::read_xlsx("data/samples_info.xlsx") %>%
  clean_names() %>%
  filter(!is.na(location_of_wes_fastq)) %>%
  select(sample_id, tumor_vs_control)

fastq_manifest <- merge(fastq_manifest, samples_manifest, by= 'sample_id') %>%
  mutate(tumor_vs_control = gsub(" ","_", tumor_vs_control))

write.table(fastq_manifest, 'outputs/GATK_preprocess_manifest.txt', quote = F, sep = "\t", row.names = F, col.names = F)
