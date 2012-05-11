\name{TBL2_MM_gen}
\alias{TBL2_MM_gen}
\title{
The target genes table for mouse miRNAs
}
\description{
This function generates target gene tables for human miRNAs base upon simple seed match toward 3' UTR of target genes}
\usage{
TBL2_MM_gen(SP="mmu",destdir="./")
}
\arguments{
  \item{SP}{Three characters string which describe mouse. Default is "mmu".}
  \item{destdir}{Directory where "id_conv" is saved. Default is "./".}
 }
\details{
This function generates target gene tables for mouse miRNAs,
  used in the MiRaGE server described in Yoshizawa et al (2011).
 It saves data frame "TBL2" as the name of "TBL2_MM" in the destination directory.
}
\references{
Yoshizawa, M., Taguchi, Y-h., and Yasuda, J. (2011), Inference of Gene Regulation via miRNAs During ES Cell
Differentiation Using MiRaGE Method. Int J Mol Sci 12(12):9265-9276
}
\examples{
\dontrun{TBL2_MM_gen()}
}

