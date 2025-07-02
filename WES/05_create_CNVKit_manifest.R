setwd(dirname(rstudioapi::getSourceEditorContext()$path))

samples_manifest <- readxl::read_xlsx("data/samples_info.xlsx") %>%
  clean_names() %>%
  filter(!is.na(location_of_wes_fastq)) %>%
  select(sample_id, tumor_vs_control)

write.table(samples_manifest, 'outputs/CNVKit_manifest.txt', quote = F, sep = "\t", row.names = F, col.names = F)
