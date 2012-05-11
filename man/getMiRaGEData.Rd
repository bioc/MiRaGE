\name{getMiRaGEData}
\alias{getMiRaGEData}
\title{
Data retrieval function for MiRaGE
}
\description{
This function retrieves target gene tables of miRNAs, gene id conversion table and miRNA consevation tables
}
\usage{
getMiRaGEData(location="local", species="MM", ID="refseq",method="mean",test="ks",
conv="conserved",species_force=T,ID_force=T,conv_force=T)
}
\arguments{
  \item{location}{Specify from where tables are retrieved. If location="local", then tables are retrieved from experimental package, "miRNATarget". If location="web", tables are downloaded from MiRaGE Server web site, "http://www.granular.com/DATA/".  Default is "local"}
  \item{species}{a string specifies target species. This should be either "MM" (mouse) or "HS" (human). Default is "MM".}
  \item{ID}{a string specifies gene ID. Default is "refseq".}
  \item{method}{a string specifies the treatment of samples. This should be one of "mean", "mixed", or "one_by_one". Default is "mean".}
  \item{test}{a string specifies statistical test for P-value computation. This should be one of "ks" (Kolmogorov-Smirnov test), "t" (t-test), or "wilcox" (Wilcoxon test). Default is "ks".}
   \item{conv}{a string specifies conservation of miRNA. This should be one of "conserved", "weak_conserve", and "all". Default is "conserved".}
   \item{species_force}{a boolean indicating whether target gene table is downloaded (TRUE) or not (FALSE). Default is TRUE.}
   \item{ID_force}{a boolean indicating whether  gene ID conversion table is downloaded (TRUE) or not (FALSE). Default is TRUE.}
   \item{conv_force}{a boolean indicating whether  miRNA conservation table is downloaded (TRUE) or not (FALSE). Default is TRUE.} 
}
\details{
This function retrieves various tables necessary  for  inference of target gene regulation by miRNAs, based upon the algorithm implemented in  the MiRaGE server described in Yoshizawa et al (2011).
}
\value{
\item{tb1}{target gene tables of miRNAs. Rownames are miRNA named and colnames are gene ids}
\item{TBL2}{Transposed matrix of tb1. This object is loaded after this function is called.}
\item{conv_id}{Conservation table of miRNAs. It is retrieved only when "conv" is not equal to "all". This object is loaded after this function is called.}
\item{id_conv}{Gene id conversion table between "refseq" ad specified gene id or probe id. This object is loaded after this function is called.} 
}
\references{
Yoshizawa, M., Taguchi, Y-h., and Yasuda, J. (2011), Inference of Gene Regulation via miRNAs During ES Cell 
Differentiation Using MiRaGE Method. Int J Mol Sci 12(12):9265-9276
}
\note{It is usually called inside the function of "MiRaGE". After the execution of getMiRaGEData, objects listed in value are automatically loaded}
\examples{
tb1 <- getMiRaGEData(species="HS")
}

