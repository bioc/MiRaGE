\name{MiRaGE}
\alias{MiRaGE}
\title{
Inference of Target Gene Regulation via MicroRNA
}
\description{
Wrapper  function computes P-values for target gene regulation via miRNAs. 
}
\usage{
MiRaGE(gene_exp,location="local",species="MM",ID="refseq",method="mean",test="ks",
conv="conserved",species_force=T,ID_force=T,conv_force=T)
}
\arguments{
  \item{gene_exp}{ExpressionSet object which stores target gene expression. featureData must include "gene_id"  and phenoData must include "sample_name". The number of columns of expression data MUST be the length of "sample_name". If not, subsetting is required.}
\item{location}{Specify from where tables are retrieved. If location="local", then tables are retrieved from experimental package, "miRNATarget". If location="web", tables are downloaded from MiRaGE Server web site, "http://www.granular.com/DATA/".  Default is "local"}
   \item{species}{a string specifies target species. This should be either "MM" (mouse) or "HS" (human).Default is "MM".}
  \item{ID}{a string specifies gene ID. Default is "refseq".}
  \item{method}{a string specifies the treatment of samples. This should be one of "mean", "mixed", or "one_by_one". Default is "mean".}
  \item{test}{a string specifies statistical test for P-value computation. This should be one of "ks" (Kolmogorov-Smirnov test), "t" (t-test), or "wilcox" (Wilcoxon test). Default is "ks".}
   \item{conv}{a string specifies conservation of miRNA. This should be one of "conserved", "weak_conserve", and "all". Default is "conserved".}
   \item{species_force}{a boolean indicating whether target gene table is downloaded (TRUE) or not (FALSE). Default is TRUE.}
   \item{ID_force}{a boolean indicating whether  gene ID conversion table is downloaded (TRUE) or not (FALSE). Default is TRUE.}
   \item{conv_force}{a boolean indicating whether  miRNA conservation table is downloaded (TRUE) or not (FALSE). Default is TRUE.} 
}
\details{
This wrapper function computes P-values for target gene regulation via miRNAs using the algorithm employed for the MiRaGE server described in Yoshizawa et al (2011). P-values are computed based upon the target gene expression change between two experimental conditions.
}
\value{
\item{P0}{data frame. The first column of "P0" contains miRNA names. From the second column to last column includes P-values. Smaller P-values mean that target genes are more likely upregulated in the later sample groups than in the former sample groups}
\item{P1}{data frame. The first column of "P1" contains miRNA names. From the second column to last column includes P-values. Smaller P-values mean that target genes are more likely downregulated in the later sample groups than  in the former sample groups}
}
\references{
Yoshizawa, M., Taguchi, Y-h., and Yasuda, J. (2011), Inference of Gene Regulation via miRNAs During ES Cell
Differentiation Using MiRaGE Method. Int J Mol Sci 12(12):9265-9276
}
\examples{
data(gene_exp)
result <- MiRaGE(gene_exp,species="HS")
result$P1[1:3,]
}

