TBL2_MM_gen <- function(SP="mmu",destdir="./")
{
#require(GenomicFeatures)
#require(BSgenome.Mmusculus.UCSC.mm9)
#require(seqinr)
txdb1 <-  makeTranscriptDbFromUCSC(genome="mm9",tablename="refGene")
UTR <- threeUTRsByTranscript(txdb1,use.name=T)
seqs <- getSeq(Mmusculus, IRanges::unlist(UTR, use.names=FALSE),as.character=TRUE)
elt <- rep(names(UTR), elementLengths(UTR))
seq2 <- sapply(split(seqs, elt), paste, collapse="")
write.fasta(sequences=lapply(as.list(seq2),s2c),file.out=paste(tempdir(),"mm9_utr3.fasta",sep="/"),names=names(seq2),nbchar=60)
download.file("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",destfile=paste(tempdir(),"mature.fa.gz",sep="/"))
cat(readLines(zz <- gzfile(paste(tempdir(),"mature.fa.gz",sep="/"))),file=paste(tempdir(),"mature.fa",sep="/"),sep="\n")
#require(Biostrings)
seqs <- read.DNAStringSet(filepath=paste(tempdir(),"mm9_utr3.fasta",sep="/"),format="fasta")
seq.ids <- names(seqs)
seq.lengths <- Biostrings::nchar(seqs)
mature <- read.RNAStringSet(filepath=paste(tempdir(),"mature.fa",sep="/"),format="fasta") 
TBL <- NULL
mature1 <-mature[grep(SP,names(mature))] #mouse data
SUB <- subseq(mature1,2,8)
for (i in c(1:length(mature1)))
{
cat(i," ")
seed.match <- DNAString(RNAString(as.character(reverseComplement(SUB[i]))))
cnt <- data.frame(ID=I(seq.ids),length=I(seq.lengths),count.7mer=I(vcountPattern(seed.match,seqs)))
TBL <- cbind(TBL,cnt[,3])
colnames(TBL)[i] <- strsplit(names(mature1[i])," ")[[1]][1]
}
TBL2<- TBL
save(file="TBL2_MM",TBL2)
}
