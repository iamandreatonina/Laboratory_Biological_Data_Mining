prova <- read.table("Complete datasets/GSE227832_RNAseq_read_counts.txt", header=T)
prova2 <- prova[332:341]
prova2[11]<-prova[1]
prova2[1]<-prova2[11]
prova2[11]<-prova[332]
colnames(prova2)[11]<-"X817_B"
colnames(prova2)[1]<-"ensembl_gene_id"

prova3<- read.table("GSE84445_Raw_counts.txt", header=T)
prova4 <- prova3[,-12:-21]

# Load the EDASeq and BioMart package
library(EDASeq)
library("biomaRt")
library(edgeR)





ensembl <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl") # create the ensembl object that points to the Hsapiens database


convert <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                 filters="hgnc_symbol", 
                 values=Controls2$Transcript_ID,
                 mart = ensembl)

# we add this info to the initial file -> use merge
complete <- merge(Controls2,convert,by.x="Transcript_ID",by.y="hgnc_symbol")

# changing the columns cause we do not need the hugo symbol
complete[1]<-complete[22]
final2<-complete[1:21]

colnames(final)[1]<-"ensembl_gene_id"

# mergin the controls -> new data begin in column 12
Controls_merg<-merge(prova2,final,by.x="ensembl_gene_id",by.y="ensembl_gene_id")

# write new table
write.csv(Controls_merg,"Complete_datasets/Control.csv", col.names = T, row.names = T, sep = ",")
