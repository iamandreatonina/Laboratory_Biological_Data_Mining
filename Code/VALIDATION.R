#install.packages("rbioapi")
#BiocManager::install("GeneOverlap")
#install.packages("epitools")
#BiocManager::install("STRINGdb")
library("GeneOverlap")
library("dplyr")
library("rbioapi")
library("epitools")
library("ggplot2")
setwd("C:/Users/lsant/OneDrive/Desktop/QCB/DataMining/project_data/rnaSeq/finals dataset")

###Retrieving IDs from OneGene expansion and STRING databse###
expansion_genes<- read.csv2("Expansion_genes.csv", header = T, sep = ",")

espansione<- as.data.frame(expansion_genes$X0)
exp_stringID <- rba_string_map_ids(ids = espansione,
                                      species = 9606)
most_connected_genes<- c("grapl","lrrc37a","llrc37a3","arl17a","stag3")
alt_gene<- c("GRAPL","LRRC37A","LRRC37A3","ARL17A","STAG3")

genes<- exp_stringID[which(exp_stringID$queryItem %in% most_connected_genes),]

genes2<- exp_stringID[which(exp_stringID$preferredName %in% alt_gene),]

genes_ID<-genes2$stringId

grapl_ID_string<- exp_stringID[which(exp_stringID$preferredName=="GRAPL"),]$stringId
llrc37a_ID_string<- exp_stringID[which(exp_stringID$preferredName=="LRRC37A"),]$stringId
llrc37a3_ID_string<- exp_stringID[which(exp_stringID$preferredName=="LRRC37A3"),]$stringId
arl17a_ID_string<- exp_stringID[which(exp_stringID$preferredName=="ARL17A"),]$stringId
stag3_ID_string<- exp_stringID[which(exp_stringID$preferredName=="STAG3"),]$stringId



###loading the expansion obtained from OneGene### #"GRAPL","LRRC37A","LRRC37A3","ARL17A","STAG3"
grapl_exp<- read.csv2("GRAPL.csv", header = T, sep = ",")
llrc37a_exp<- read.csv2("LRRC37A.csv", header = T, sep = ",")
llrc37a3_exp<- read.csv2("LRRC37A3.csv", header = T, sep = ",")
arl17a_exp<- read.csv2("ARL17A.csv", header = T, sep = ",")
stag3_exp<- read.csv2("STAG3.csv", header = T, sep = ",")



###creating the STRING expansions###

#grapl
grapl_exp_string<-rba_string_interactions_network(
  ids = grapl_ID_string,
  species =  9606,
  required_score = NULL,
  network_type = "functional",

  use_query_labels = FALSE
)

#llrc37a
llrc37a_exp_string<-rba_string_interactions_network(
  ids = llrc37a_ID_string,
  species =  9606,
  required_score = NULL,
  network_type = "functional",
  use_query_labels = FALSE
)

#llrc37a3
llrc37a3_exp_string<-rba_string_interactions_network(
  ids = llrc37a3_ID_string,
  species =  9606,
  required_score = NULL,
  network_type = "functional",
  use_query_labels = FALSE
)

#arl17a
arl17a_exp_string<-rba_string_interactions_network(
  ids = arl17a_ID_string,
  species =  9606,
  required_score = NULL,
  network_type = "functional",
  use_query_labels = FALSE
)

#stag3
stag3_exp_string<-rba_string_interactions_network(
  ids = stag3_ID_string,
  species =  9606,
  required_score = NULL,
  network_type = "functional",
  use_query_labels = FALSE
)

###overlapping the Expansions and testing them###

#grapl
grapl_exp <- dplyr::mutate_all(grapl_exp, .funs = toupper)
count<- grapl_exp_string %>% dplyr::filter(grapl_exp_string$preferredName_A %in% grapl_exp$X0)

overlap_grapl<-newGeneOverlap(grapl_exp$X0, grapl_exp_string$preferredName_A)
go_grapl.obj<-testGeneOverlap(overlap_grapl)
getContbl(go_grapl.obj)
print (go_grapl.obj)


#llrc37a
llrc37a_exp <- dplyr::mutate_all(llrc37a_exp, .funs = toupper)
count<- llrc37a_exp_string %>% dplyr::filter(llrc37a_exp_string$preferredName_A %in% llrc37a_exp$X0)

overlap_llrc37a<-newGeneOverlap(llrc37a_exp$X0, llrc37a_exp_string$preferredName_A)
go_llrc37a.obj<-testGeneOverlap(overlap_llrc37a)
getContbl(go_llrc37a.obj)
print (go_llrc37a.obj)




#llrc37a3
llrc37a3_exp <- dplyr::mutate_all(llrc37a3_exp, .funs = toupper)
count<- llrc37a3_exp_string %>% dplyr::filter(llrc37a3_exp_string$preferredName_A %in% llrc37a3_exp$X0)

overlap_llrc37a3<-newGeneOverlap(llrc37a3_exp$X0, llrc37a3_exp_string$preferredName_A)
go_llrc37a3.obj<-testGeneOverlap(overlap_llrc37a3)
getContbl(go_llrc37a3.obj)
print (go_llrc37a3.obj)

#arl17a
arl17a_exp <- dplyr::mutate_all(arl17a_exp, .funs = toupper)
count<- arl17a_exp_string %>% dplyr::filter(arl17a_exp_string$preferredName_A %in% arl17a_exp$X0)

overlap_arl17a<-newGeneOverlap(arl17a_exp$X0, arl17a_exp_string$preferredName_A)
go_arl17a.obj<-testGeneOverlap(overlap_arl17a)
getContbl(go_arl17a.obj)
print (go_arl17a.obj)

#stag3
stag3_exp <- dplyr::mutate_all(stag3_exp, .funs = toupper)
count<- stag3_exp_string %>% dplyr::filter(stag3_exp_string$preferredName_A %in% stag3_exp$X0)

overlap_stag3<-newGeneOverlap(stag3_exp$X0, stag3_exp_string$preferredName_A)
go_stag3.obj<-testGeneOverlap(overlap_stag3)
getContbl(go_stag3.obj)
print (go_stag3.obj)

###GRAPL vs all expansion###
grapl_exp <- dplyr::mutate_all(grapl_exp, .funs = toupper)
espansione <- dplyr::mutate_all(espansione, .funs = toupper)
espansione_l <- as.list(espansione$`expansion_genes$X0`)
count<- grapl_exp_string %>% dplyr::filter(grap_exp_string$preferredName_A %in% espansione_l)

overlap_grapl_all<-newGeneOverlap(espansione_l, grapl_exp_string$preferredName_A)
go_grapl_all.obj<-testGeneOverlap(overlap_grapl_all)
getContbl(go_grapl_all.obj)
print (go_grapl_all.obj)

###grap as alternative to grapl
grap_ID_string<- exp_stringID[which(exp_stringID$preferredName=="GRAP"),]$stringId
grap_exp_string<-rba_string_interactions_network(
  ids = grap_ID_string,
  species =  9606,
  required_score = NULL,
  network_type = "functional",
  use_query_labels = FALSE
)

grapl_exp <- dplyr::mutate_all(grapl_exp, .funs = toupper)
count<- grap_exp_string %>% dplyr::filter(grap_exp_string$preferredName_A %in% grapl_exp$X0)

overlap_grap<-newGeneOverlap(grapl_exp$X0, grap_exp_string$preferredName_A)
go_grap.obj<-testGeneOverlap(overlap_grap)
getContbl(go_grap.obj)
print (go_grap.obj)


