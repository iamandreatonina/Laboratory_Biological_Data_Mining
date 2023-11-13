library(tidyverse)

file_2018 <- read.table('all_phase2_target_2018_pub/data_mrna_seq_rpkm.txt',header = T,sep = '\t')
genes_file <- readxl::read_xlsx('Human-specific.xlsx')

file_clinical_status_sample<- read.table('all_phase2_target_2018_pub/data_clinical_sample.txt',header = T,sep = '\t')
file_clinical_status_patient <- read.table('all_phase2_target_2018_pub/data_clinical_patient (copy).txt', header = T ,sep = '\t',fill = T)


file_filtered <- file_2018 %>% filter(Hugo_Symbol %in% genes_file$`Gene Name`)


