\name{MM_conv_id}
\alias{MM_conv_id}
\title{
miRNA conservation table generation function
}
\description{
This function generates conv_id, which describes conservation information of each mouse miRNA, from the scratch}
\usage{
MM_conv_id(taxid=10090,species="mmu",destdir="./")
}
\arguments{
  \item{taxid}{taxon id for human. Default is 11090}
  \item{species}{Three characters string which describe human. Default is "mmu".}
  \item{destdir}{Directory where "MM_conv_id" is saved. Default is "./".}
 }
\details{
This function computes conv_id, which describes conservation information of each mouse miRNA,
  used in the MiRaGE server described in Yoshizawa et al (2011).
 It saves data frame conv_id as the name of "MM_conv_id" in the destination directory.
}
\references{
Yoshizawa, M., Taguchi, Y-h., and Yasuda, J. (2011), Inference of Gene Regulation via miRNAs During ES Cell
Differentiation Using MiRaGE Method. Int J Mol Sci 12(12):9265-9276
}
\examples{
\dontrun{MM_conv_id()}
}

