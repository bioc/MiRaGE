getMiRaGEData <- function(location="local", species="MM", ID="refseq",method="mean",test="ks",conv="conserved",species_force=T,ID_force=T,conv_force=T)
{
if (location=="local")
 {
 if(!require(miRNATarget, quietly=TRUE)) 
  {
  cat("Package:miRNATarget not available attempting to biocLite\n")
  source("http://bioconductor.org/biocLite.R")
  biocLite("miRNATarget")
  require(miRNATarget)
  }
 } else if (location=="web") {
  loadurl <- function(x)
  {
  con <- url(paste("http://granular.com/MiRaGE/DATA/",x,sep=""))
  load(con,envir = parent.frame(n=3))
  close(con)
  }
 }



if (species_force)
{
if (species=="MM")
{

if (location=="local")
{
data(TBL2_MM)
}else if (location=="web"){
loadurl("TBL2_MM")
}
tb1 <- t(TBL2)
} else if (species=="HS") {

if (location=="local")
{
data(TBL2_HS)
}else if (location=="web"){
loadurl("TBL2_HS")
}
tb1 <- t(TBL2)
} else {
print("Wrong species!")
return()
}
} else{
if ("TBL2" %in% objects(envir = parent.frame(n=1))) {
tb1 <- t(TBL2)
} #else {
#print("No target gene table was loaded. Specify species")
#return()
#}
}
if (ID_force)
{
if (ID != "refseq")
{

if (location=="local")
{
data(list=paste(species,"_refseq_to_",ID,sep=""))
} else if (location=="web"){
loadurl(paste("./id_conv/",species,"/",species,"_refseq_to_",ID,sep=""))
}

index <- id_conv[match(colnames(tb1),id_conv[,1]),2]
index2 <- is.na(index) | index==""
tb1 <- tb1[,!index2]
colnames(tb1) <- index[!index2]
}
}else{
if ("id_conv" %in% objects(envir = parent.frame(n=1))) {
index <- id_conv[match(colnames(tb1),id_conv[,1]),2]
index2 <- is.na(index) | index==""
tb1 <- tb1[,!index2]
colnames(tb1) <- index[!index2]
} #else {
#print("No  gene ID conservation table was downloaded. Specify ID")
#return()
#}
}


if (conv_force)
{
if (conv!="all")
{

if (location=="local")
{
data(list=paste(species,"_conv_id",sep=""))
} else if (location=="web"){
loadurl(paste(species,"_conv_id",sep=""))
}

if (conv=="conserved") {
index <- rownames(tb1) %in% conv_id[ conv_id[,2]==2,1]
} else if (conv=="weak_conserv") {
index <- rownames(tb1) %in% conv_id[ conv_id[,2]!=0,1]
}
tb1 <- tb1[index,]
}
}else{
if ("conv_id" %in% objects(envir = parent.frame(n=1))) {
if (conv=="conserved") {
index <- rownames(tb1) %in% conv_id[ conv_id[,2]==2,1]
} else if (conv=="weak_conserv") {
index <- rownames(tb1) %in% conv_id[ conv_id[,2]!=0,1]
}
tb1 <- tb1[index,]
} #else {
#print("No  gene ID conservation table was downloaded. Specify conservation")
#return()
#}
}

return(tb1)
}
