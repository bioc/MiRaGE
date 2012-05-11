\name{TBL2_HS_gen}
\alias{TBL2_HS_gen}
\title{
The target genes table for human miRNAs
}
\description{
This function generates target gene tables for human miRNAs base upon simple seed match toward 3' UTR of target genes}
\usage{
TBL2_HS_gen(SP="hsa",destdir="./")
}
\arguments{
  \item{SP}{Three characters string which describe human. Default is "hsa".}
  \item{destdir}{Directory where "TBL2_HS" is saved. Default is "./".}
 }
\details{
This function generates target gene tables for human miRNAs,
  used in the MiRaGE server described in Yoshizawa et al (2011).
 It saves data frame "TBL2" as the name of "TBL2_HS" in the destination directory.
}
\references{
Yoshizawa, M., Taguchi, Y-h., and Yasuda, J. (2011), Inference of Gene Regulation via miRNAs During ES Cell
Differentiation Using MiRaGE Method. Int J Mol Sci 12(12):9265-9276
}
\examples{
\dontrun{TBL2_HS_gen()}
}

