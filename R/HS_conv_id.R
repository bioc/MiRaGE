HS_conv_id <- function(taxid=9606,species="hsa",destdir="./")
{
download.file("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",destfile=paste(tempdir(),"mature.fa.gz",sep="/"))
cat(readLines(zz <- gzfile(paste(tempdir(),"mature.fa.gz",sep="/"))),file=paste(tempdir(),"mature.fa",sep="/"),sep="\n")
download.file("http://www.targetscan.org/vert_61/vert_61_data_download/miR_Family_Info.txt.zip",destfile=paste(tempdir(),"miR_Family_Info.txt.zip",sep="/"))
cat(readLines(zz <- unz(paste(tempdir(),"miR_Family_Info.txt.zip",sep="/"),filename="miR_Family_Info.txt")),file=paste(tempdir(),"miR_Family_Info.txt",sep="/"),sep="\n")
#require(seqinr)
x <- read.fasta(paste(tempdir(),"mature.fa",sep="/"))
x1 <- read.csv(paste(tempdir(),"miR_Family_Info.txt",sep="/"),sep="\t")
miRNA<- names(x)[grep(species,names(x))]
Annot <- lapply(x,attr,"Annot")
Accession <- unlist(Annot)[grep(species,names(x))]
Accession <- lapply(strsplit(Accession," "),"[",2)
index <- match(Accession,x1[x1[,3]==taxid,7])
conv_id <- data.frame(miRNA,x1[x1[,3]==taxid,6][index])
conv_id[is.na(conv_id[,2]),2] <-0
save(file=paste(destdir,"HS_conv_id",sep=""),conv_id)
}
