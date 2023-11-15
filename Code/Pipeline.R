# BiocManager::install('sva')
library(sva)
library(tibble)

##### Batch effect correction 

#Load datasets 
Tumor <- read.csv('Tumor_dataframe.csv',sep =',',header = T)
Control <- read.csv('Controls.csv',sep = ',',header = T)

Tumor_2 <- as.matrix(sapply(Tumor[2:641], as.numeric))
Control_2 <- as.matrix(sapply(Control[2:31], as.numeric))

#creation of batch for tumor and control, so creation of the vectors for batch separation 
batch_tumor <- c(rep(1,38),rep(2,173),rep(3,321),rep(4,108))
batch_control <- c(rep(1,10),rep(2,20))

# application of Combat-Seq and creation of adjusted dataframes 
tumor_adjusted <- as.data.frame(ComBat_seq(Tumor_2,batch = batch_tumor,group = NULL))
control_adjusted <- as.data.frame(ComBat_seq(Control_2, batch = batch_control, group = NULL))

# adding the ensembl_gene_id column 
control_adjusted <- add_column(control_adjusted,'ensembl_gene_id' =Control$ensembl_gene_id, .before = 'X817_T')
tumor_adjusted <- add_column(tumor_adjusted, 'ensembl_gene_id' = Tumor$ensembl_gene_id, .before = 'ALL33')

##### Normalization with edgeR package 

# We use TMM method , which is a normalization method intra and inter-sample and we create CPM matrices 
library(edgeR)
install.packages('DESeq2')
library(DESeq2)

# We don't have a precise threshold to eliminate the low expressed genes, so we know that DESeq2 set a 
#threshold based on the expression of our data by doing result()

# set the dataframe more easier for us to use 
control_adjusted1 <- control_adjusted %>% column_to_rownames('ensembl_gene_id')
tumor_adjusted1 <- tumor_adjusted %>% column_to_rownames('ensembl_gene_id')

# Presence of a duplicate gene in the Tumor dataset, think what to do 
Tumor$ensembl_gene_id[duplicated(Tumor$ensembl_gene_id)]
duplicate <- tumor_adjusted %>% dplyr::filter(tumor_adjusted$ensembl_gene_id == 'ENSG00000230417')

# first we find the DGEList object 

edge_c_control <- DGEList(counts = control_adjusted1)
edge_c_tumor <- DGEList(counts = tumor_adjusted1)






