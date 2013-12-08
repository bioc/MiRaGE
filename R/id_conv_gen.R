id_conv_gen <- function(SP="MM",gene_id="all",destdir="./")
{
#require("biomaRt")
ensembl = useMart("ensembl")
listDatasets(ensembl)
if (SP=="HS") {ensembl <-  useDataset("hsapiens_gene_ensembl",mart=ensembl)}
if (SP=="MM") {ensembl <-  useDataset("mmusculus_gene_ensembl",mart=ensembl)}
features <- data.frame(listAttributes(ensembl, page = "feature_page"))
write.table(file=paste(tempdir(),"/features.csv",sep=""),features,quote=F,row.names=F,sep="\t")
features <- read.csv(paste(tempdir(),"/features.csv",sep=""),sep="\t")
counts<-NULL
if (gene_id=="all") {
i1 <- 1
i2 <- dim(features)[1]
} else {
i1 <- grep(gene_id,features[,1])
i2 <- i1
}
for (i in c(i1:i2)) 
#for (i in c(1:2)) 
{
cat(i," ")
id <- features[i,1]
if (id == "refseq_mrna" | id == "refseq_ncrna") next
id_conv <- getBM(attributes = c("refseq_mrna",as.character(id)),mart = ensembl)
id_conv2 <- getBM(attributes = c("refseq_ncrna",as.character(id)),mart = ensembl)
colnames(id_conv)[1] <- "refseq_mrna_ncrna"
colnames(id_conv2)[1] <- "refseq_mrna_ncrna"
id_conv <- rbind(id_conv,id_conv2)
index <- id_conv[,1]!="" & id_conv[,2]!=""
counts <- c(counts,sum(index,na.rm=T))
id_conv<-id_conv[index,]
save(file=paste(destdir,SP,"_refseq_to_",id,sep=""),id_conv)
}
}
