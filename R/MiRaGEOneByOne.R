MiRaGEOneByOne <- function(x_gene,tb1,gp,s_test_less,s_test_greater)
{
for ( i0 in c(1:dim(gp)[1]))
{
if (i0 < dim(gp)[1])
{
for ( i1 in c((i0+1):dim(gp)[1]))
{
if (i1 <= dim(gp)[1])
{

c0 <- c(1:dim(gp)[2])[gp[i0,]]+1
c1 <- c(1:dim(gp)[2])[gp[i1,]]+1
P0 <- NULL
P1 <- NULL
for ( i  in c0)
{
for ( j in c1)
{
x11 <- log(x_gene[,i]/x_gene[,j])
P00 <- NULL
P11 <- NULL

for (k in c(1:dim(tb1)[1]))
{
P00<-  c(P00,s_test_less(x11[match(colnames(tb1)[tb1[k,]>0],x_gene[,1])],x11[match(colnames(tb1)[tb1[k,]==0],x_gene[,1])]))
P11<-  c(P11,s_test_greater(x11[match(colnames(tb1)[tb1[k,]>0],x_gene[,1])],x11[match(colnames(tb1)[tb1[k,]==0],x_gene[,1])]))
}
P0 <- cbind(P0,P00)
P1 <- cbind(P1,P11)
}
}

P0 <- data.frame(rownames(tb1),P0)
P1 <- data.frame(rownames(tb1),P1)

names<-NULL
for ( i  in c0) 
{
for ( j in c1)
{
names <- c(names,paste(colnames(x_gene)[i],"vs",colnames(x_gene)[j]))
}
}
colnames(P0) <- c("Refseq",names)
colnames(P1) <- c("Refseq",names)
}}
}}
#save(file=paste("/tmp/",path,"/P0",sep=""),P0)
#save(file=paste("/tmp/",path,"/P1",sep=""),P1)
#source("results2.R")
return(list(P0=P0,P1=P1))
}

