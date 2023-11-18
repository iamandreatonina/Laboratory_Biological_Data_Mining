# BiocManager::install('sva')
library(sva)
library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(stringr)

setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Datasets_finals")

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

# Save the CPM table 
write.csv(CPM_tumor,file = 'CPM_Tumor_dataframe.csv',row.names = T)
write.csv(CPM_control,file = 'CPM_Control_dataframe.csv',row.names = T)

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

###### Differential gene expression analysis 

total_adjusted <- merge(control_adjusted,tumor_adjusted,by='ensembl_gene_id')
total_adjusted1 <- total_adjusted %>% column_to_rownames('ensembl_gene_id')

# creating a dataframe containing the info on the samples, this is needed to be able to perform the DGE, we set the conditions of the samples as healty (H) and tumoral (T)
info_sample_1<-data.frame("sample"=colnames(total_adjusted1))
rownames(info_sample_1)<-info_sample_1$sample
info_sample_2<-as.data.frame(str_split(string=info_sample_1$sample, pattern="_", simplify=T))
colnames(info_sample_2)<-c("condition","replicate")
info_samples<-cbind(info_sample_1, info_sample_2[1:2])
info_samples$condition<-c(rep("H",30),rep("T",640))
info_samples$replicate<-c(rep(1,670))


# Now we can create the DGEList object
edge_c_total <- DGEList(counts = total_adjusted1, group=info_samples$condition, samples=info_samples, genes=total_adjusted1)
edge_n_total <- calcNormFactors(edge_c_total,method = 'TMM')

# We create the cpm table
cpm_table <-as.data.frame(round(cpm(edge_n_total),2)) # the library size is scaled by the normalization factor

# Here we define the experimental design matrix, we build a model with no intercept also we have two varaibles, one for each condition 
# 1 for control and 2 for tumor 
design <- model.matrix(~0+group, data = edge_n_total$samples)
colnames(design) <- levels(edge_n_total$samples$group)
rownames(design) <- edge_n_total$samples$sample

# Calculate dispersion and fit the result with edgeR (necessary for differential expression analysis)
edge_d_total <- estimateDisp(edge_n_total,design)

# Fit the data we model the data using a negative binomial distribution
edge_f<-glmQLFit(edge_d_total, design)

# Definition of the contrast (conditions to be compared)
contro <- makeContrasts("H-T", levels=design)

# Fit the model with generalized linear models
edge_t <- glmQLFTest(edge_f,contrast=contro)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change
DEGs <- as.data.frame(topTags(edge_t,n=18203,p.value = 0.01,sort.by = "logFC"))

# We add a new column to the DEGs dataframe called class.
# Used to express the values of the fold change of the transcripts.
# The selection is based on the log fold change ratio (>1.5 for up-regulated genes and < (-1.5) for down-regulated genes)
# and a log CPM (>1 for both cases). From the contingency table of our DEGs we can see that the up regulated genes
# correspond to the 3.7% of the total and the down regulated are the 16% of the total.

DEGs$class <- '='
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC > 1.5)] = '+'
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC < (-1.5))] = '-'
DEGs <- DEGs[order(DEGs$logFC, decreasing = T),] # we order based on the fold change

table(DEGs$class)
# -:2942,+:683,=:7669


# Let`s check how many human specific genes we have in the up regulated and down regulated genes
#  We have 110 down-reg HS genes and 18 up-regulated HS genes
DEGs_Hsgenes <- DEGs %>% dplyr::filter(rownames(DEGs) %in% Human_genes$`Ensembl ID`)
Up_HSgenes <- DEGs[DEGs$class=='+',] %>% dplyr::filter(rownames(DEGs[DEGs$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HSgenes <- DEGs[DEGs$class=='-',] %>% dplyr::filter(rownames(DEGs[DEGs$class=='-',]) %in% Human_genes$`Ensembl ID`) 


# Display the results using a volcano plot (x-axes: log FoldChange, y-axes: inverse function of the p-value).
# We can see the most significant DEGs colored in green, which are genes that surpass a threshold set on both the p-value
# and the Fold Change.
jpeg(filename = '../images/Vulcano_plot_DEGs.jpeg')
input_df<-DEGs
xlabel<- "log2 FC control vs case"
ylabel<-"-log10 p-value"
par(fig=c(0,1, 0,1), mar=c(4,4,1,2), mgp=c(2, 0.75,0))
plot(DEGs$logFC,-log(DEGs$PValue, base=10), xlab=xlabel,ylab = ylabel, col=ifelse(DEGs$class=="=", "grey70", "olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot") %>% 
abline(v = 0, lty = 2, col="grey20")
dev.off()

######### vulcano hs
jpeg(filename = '../images/Vulcano_plot_DEGsHS.jpeg')
input_df<-DEGs_Hsgenes
xlabel<- "log2 FC control vs case"
ylabel<-"-log10 p-value"
par(fig=c(0,1, 0,1), mar=c(4,4,1,2), mgp=c(2, 0.75,0))
plot(DEGs_Hsgenes$logFC,-log(DEGs_Hsgenes$PValue, base=10), xlab=xlabel,ylab = ylabel, col=ifelse(DEGs_Hsgenes$class=="=", "grey70", "olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot") %>% 
  abline(v = 0, lty = 2, col="grey20")
dev.off()

# We can also represent the genes using a heatmap.
# A clustering process is operated. We plot only up or down expressed genes using data from both the normalized CPM and
# the log transformation of the CPM table.

col <- rep('chartreuse4', 670)
col[which(info_samples$condition == 'T')] <- 'burlywood3' 
pal <- c('blue','white','red')
pal <- colorRampPalette(pal)(670)
DEGs_selected <- DEGs %>% dplyr::filter(DEGs$class != '=')
jpeg(filename = '../images/Heatmap_plot_DEGs.jpeg')
heatmap(as.matrix(cpm_table[which(rownames(cpm_table) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()

#for improving the clusterization we set the cpm table as logarithmic
cpm_table_log <- as.data.frame(round(log10(cpm(edge_n_total)+1),2))
jpeg(filename = '../images/Heatmap_plot_DEGs_log.jpeg')
heatmap(as.matrix(cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()

############ heatmap human specific

col <- rep('chartreuse4', 670)
col[which(info_samples$condition == 'T')] <- 'burlywood3' 
pal <- c('blue','white','red')
pal <- colorRampPalette(pal)(670)
DEGs_selected <- DEGs_Hsgenes %>% dplyr::filter(DEGs_Hsgenes$class != '=')
jpeg(filename = '../images/Heatmap_plot_DEGsHS.jpeg')
heatmap(as.matrix(cpm_table[which(rownames(cpm_table) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()

#for improving the clusterization we set the cpm table as logarithmic
cpm_table_log <- as.data.frame(round(log10(cpm(edge_n_total)+1),2))
jpeg(filename = '../images/Heatmap_plot_DEGs_logHS.jpeg')
heatmap(as.matrix(cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()


##### PCA analysis 
PCA_cpm_log <- cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),]

# # we need to have both for the columns and the row a variance different from zero (because divide for the varaince )
PCA_cpm_log_filtered<-PCA_cpm_log[,which(apply(PCA_cpm_log, 2, var) != 0)]
PCA_cpm_log_filtered<- PCA_cpm_log_filtered[which(apply(PCA_cpm_log_filtered, 1, var) != 0),]
color<- c(rep('darkgreen',30),rep('indianred',640))
# # PCA plot
data.PC <- prcomp(t(PCA_cpm_log_filtered),scale. = T)
jpeg(filename = '../images/PCA_plot_DEGs_log_HS.jpeg')
plot(data.PC$x[,1:2],col=color,pch = 19) # presence of an outlier for PC1 over 500000
dev.off()

# # PCA plot of tumor only
data.PC.tumor <- prcomp(t(PCA_cpm_log_filtered[31:670]),scale. = T)
jpeg(filename = '../images/PCA_plot_DEGs_log_tumor.jpeg')
plot(data.PC.tumor$x[,1:2],pch = 19) # presence of an outlier for PC1 over 500000
dev.off()

##### Partitioning around medoids, need to also to install cmake
# install.packages('factoextra')
# install.packages('cluster')
library(factoextra)
library(cluster)

# for control-tumor -> 2 clusters as seen from graphs under
fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10)
fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10,method = 'gap_stat')+ theme_classic()
fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10, method = "wss")

# For subtypes of tumors
fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15)
fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15,method = 'gap_stat')+ theme_classic()
fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15, method = "wss")

# PAM on control-tumor
pam1<-eclust(data.PC$x, "pam", k=9)

# PAM tumors subypes
pam2<-eclust(data.PC.tumor$x, "pam", k=8)

##### hierarchical clustering

#calculate distances between observations and create a simple dendogram
dm <- dist(data.PC$x)
hc <- hclust(dm,method = 'average')
plot(hc,hang =-1)
rect.hclust(hc,k=2,border = 'red')
clust.vec.2<-cutree(hc,k=2)
fviz_cluster(list(data=data.PC$x, cluster=clust.vec.2))

# same for tumor subtypes
dm2 <- dist(data.PC.tumor$x)
hc2 <- hclust(dm2,method = 'average')
plot(hc2,hang =-1)
rect.hclust(hc2,k=2,border = 'red')
clust.vec<-cutree(hc2,k=8)
fviz_cluster(list(data=data.PC.tumor$x, cluster=clust.vec))

clusters <- mutate(cpm_table_log[31:670],cluster =clust.vec)
