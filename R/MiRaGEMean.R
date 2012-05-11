MiRaGEMean <- function(x_gene,tb1,gp,s_test_less,s_test_greater)
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

xx11 <- NULL
for (i in c0){
for (j in c1){
xx11 <- cbind(xx11,log(x_gene[,i]/x_gene[,j]))
}}

x11 <- rowMeans(xx11)

P0 <- NULL
P1 <- NULL
for (k in c(1:dim(tb1)[1]))
{
xa <- x11[match(colnames(tb1)[tb1[k,]>0],x_gene[,1])]
xb <- x11[match(colnames(tb1)[tb1[k,]==0],x_gene[,1])]
if (length(xa)>1 & length(xb)>1) 
{
P0<-c(P0,s_test_less(xa,xb))
P1<-c(P1,s_test_greater(xa,xb))
}
}

P0 <- data.frame(rownames(tb1),P0)
P1 <- data.frame(rownames(tb1),P1)
colnames(P0) <- c("Refseq","mixed")
colnames(P1) <- c("Refseq","mixed")
return(list(P0=P0,P1=P1))
}}
}}

   
}
