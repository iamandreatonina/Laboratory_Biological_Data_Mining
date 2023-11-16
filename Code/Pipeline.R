# BiocManager::install('sva')
library(sva)
library(tibble)
library(dplyr)
library(readxl)
library(ggplot2)


##### Upload huamn specific genes 

Human_genes <- readxl::read_xlsx('Human-specific.xlsx')

##### Batch effect correction 

#Load datasets 
Tumor <- read.csv('Tumor_dataframe.csv',sep =',',header = T)
Control <- read.csv('Controls.csv',sep = ',',header = T)

# We found out a duplicated ensembl_gene_id due to the fact there isn't a 1 to 1 mapping from ensembl_gene_id and hugo_symbols
# so we are going to eliminate the less informative one 

duplicato <- Tumor$ensembl_gene_id[duplicated(Tumor$ensembl_gene_id)]
# sum <- Tumor %>% dplyr::filter(Tumor$ensembl_gene_id == duplicato) 
# rowSums(sum[2:641]) # the first one is the most informative so we use distinct()
Tumor <- distinct(Tumor,ensembl_gene_id,.keep_all =T )

Tumor_2 <- as.matrix(sapply(Tumor[2:641], as.numeric))
Control_2 <- as.matrix(sapply(Control[2:31], as.numeric))

PreCombat_control_df <- tidyr::gather(as.data.frame(Control_2),key = 'sample',value = 'read_number')
PreCombat_tumor_df_subset <- tidyr::gather(as.data.frame(Tumor_2)[1:20],key = 'sample',value = 'read_number')

jpeg(filename = '../images/control_Pre_Combat_boxplot.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkgreen',fill = 'olivedrab',alpha = 0.7, mapping = aes(sample,read_number+1),data = PreCombat_control_df)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

# in this case there are too many samples so we are gonging to plot 20 samples instead of 640 
jpeg(filename = '../images/tumor_Pre_Combat_boxplot_subset.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkred',fill = 'indianred',alpha = 0.5, mapping = aes(sample,read_number+1),data = PreCombat_tumor_df_subset, width = 0.5)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

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
# install.packages('DESeq2')
library(DESeq2)
library(tidyverse)

# Let`s check how many human specific genes we have in our dataset
HSgenes_tumor <- tumor_adjusted %>% dplyr::filter(tumor_adjusted$ensembl_gene_id %in% Human_genes$`Ensembl ID`) 
# the result is that of 873 human specific genes there are present 603 in tumor  

HSgenes_control <- control_adjusted %>% dplyr::filter(control_adjusted$ensembl_gene_id %in% Human_genes$`Ensembl ID`) 
# the result is that of 873 human specific genes there are present 498 in control 


# We don't have a precise threshold to eliminate the low expressed genes, so we know that DESeq2 set a 
#threshold based on the expression of our data by doing result()

# set the dataframe more easier for us to use 
control_adjusted1 <- control_adjusted %>% column_to_rownames('ensembl_gene_id')
tumor_adjusted1 <- tumor_adjusted %>% column_to_rownames('ensembl_gene_id')

#let's have a look at the data using a boxplot, we will make a comparison after the normalization
Pre_control_df <- tidyr::gather(control_adjusted1,key = 'sample',value = 'read_number')
Pre_tumor_df_subset <- tidyr::gather(tumor_adjusted1[1:20],key = 'sample',value = 'read_number')

jpeg(filename = '../images/control_Pre_TMM_boxplot.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkgreen',fill = 'olivedrab',alpha = 0.7, mapping = aes(sample,read_number+1),data = Pre_control_df)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

# in this case there are too many samples so we are gonging to plot 20 samples instead of 640 
jpeg(filename = '../images/tumor_Pre_TMM_boxplot_subset.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkred',fill = 'indianred',alpha = 0.5, mapping = aes(sample,read_number+1),data = Pre_tumor_df_subset, width = 0.5)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

# first we find the DGEList object 

edge_c_control <- DGEList(counts = control_adjusted1)
edge_c_tumor <- DGEList(counts = tumor_adjusted1)

# normalize with the edgeR package using the TMM method, which apply and inter and intra normalization of the data, both for the controls and the tumor 
edge_n_control <- calcNormFactors(edge_c_control,method = 'TMM') 
edge_n_tumor <- calcNormFactors(edge_c_tumor,method = 'TMM')

# from that we create a CPM table with normalized expression values 
CPM_control <- as.data.frame(round(cpm(edge_n_control),2)) 
CPM_tumor <-  as.data.frame(round(cpm(edge_n_tumor),2))

CPM_control_df <- tidyr::gather(CPM_control,key = 'sample',value = 'CPM')
CPM_tumor_df <- tidyr::gather(CPM_tumor,key = 'sample',value = 'CPM')
CPM_tumor_df_toplot <- tidyr::gather(CPM_tumor[1:20],key = 'sample',value = 'CPM')

jpeg(filename = '../images/control_TMM_boxplot.jpeg')
ggplot() +
         geom_boxplot(colour = 'darkgreen',fill = 'olivedrab',alpha = 0.7, mapping = aes(sample,CPM+1),data = CPM_control_df)+
         theme_bw()+
         scale_x_discrete(guide = guide_axis(angle = 90))+
         scale_y_log10()
dev.off()

# in this case there are too many samples so we are gonging to plot 20 samples instead of 640 
jpeg(filename = '../images/tumor_TMM_boxplot_subset.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkred',fill = 'indianred',alpha = 0.5, mapping = aes(sample,CPM+1),data = CPM_tumor_df_toplot, width = 0.5)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

##### PCA analysis 

CPM_table_total <- as.data.frame(merge(CPM_control,CPM_tumor,by = 0))
color<- c(rep('darkgreen',29),rep('indianred',640))
# we need to have both for the columns and the row a variance different from zero (because divide for the varaince )
CPM_table_total_filtered<-CPM_table_total[ , which(apply(CPM_table_total, 2, var) != 0)] 
CPM_table_total_filtered<- CPM_table_total_filtered[which(apply(CPM_table_total_filtered, 1, var) != 0),]
CPM_table_total_filtered_PCA<-CPM_table_total_filtered[2:670]
# test <- CPM_table_total_filtered_PCA +1

data.PC <- prcomp(t(test),scale. = F)
plot(data.PC$x[,1:2],col=color,pch = 19,outlie)
