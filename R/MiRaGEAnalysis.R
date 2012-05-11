MiRaGEAnalysis <- function(gene_exp, tb1, method="mean", test="ks")
{

x_gene <- data.frame(fData(gene_exp)$gene_id,exprs(gene_exp))
colnames(x_gene) <- c("gene",pData(gene_exp)$sample_name)

index <- x_gene[,1] %in% colnames(tb1)
x_gene <- x_gene[index,]
index <-colnames(tb1) %in% x_gene[,1]
tb1<- tb1[,index]
#save(file=paste("/tmp/",path,"/tb1",sep=""),tb1)
#save(file=paste("/tmp/",path,"/x_gene",sep=""),x_gene)

cl <-colnames(x_gene)
df <-t(data.frame(strsplit(cl[2:length(cl)],".",fixed=T)))
sname <- unique(df[,1])
gp <- outer(sname,df[,1],"==")

if (test=="t")
{
s_test_less <- function(x1,x2){t.test(x1,x2,alternative="less")$p.value}
s_test_greater <- function(x1,x2){t.test(x1,x2,alternative="greater")$p.value}
}else if (test=="wilcox"){
s_test_less <- function(x1,x2){wilcox.test(x1,x2,alternative="less")$p.value}
s_test_greater <- function(x1,x2){wilcox.test(x1,x2,alternative="greater")$p.value}
}else if (test=="ks")
{
s_test_less <- function(x1,x2){ks.test(x1,x2,alternative="greater")$p.value}
s_test_greater <- function(x1,x2){ks.test(x1,x2,alternative="less")$p.value}
}
  if (method == "mean")
     return(MiRaGEMean(x_gene,tb1,gp,s_test_less,s_test_greater))
    else if (method == "mixed")
      return(MiRaGEMixed(x_gene,tb1,gp,s_test_less,s_test_greater))
    else if (method == "one_by_one")
       return(MiRaGEOneByOne(x_gene,tb1,gp,s_test_less,s_test_greater))
}
