\name{HS_conv_id}
\alias{HS_conv_id}
\title{
miRNA conservation table generation  function
}
\description{
This function generates conv_id, which describes conservation information of each human miRNA, from the scratch}
\usage{
HS_conv_id(taxid=9606,species="hsa",destdir="./")
}
\arguments{
  \item{taxid}{taxon id for human. Default is 9606}
  \item{species}{Three characters string which describe human. Default is "hsa".}
  \item{destdir}{Directory where "HS_conv_id" is saved. Default is "./".}
 }
\details{
This function computes conv_id, which describes conservation information of each human miRNA,
  used in the MiRaGE server described in Yoshizawa et al (2011).
 It saves data frame conv_id as the name of "HS_conv_id" in the destination directory.
}
\references{
Yoshizawa, M., Taguchi, Y-h., and Yasuda, J. (2011), Inference of Gene Regulation via miRNAs During ES Cell
Differentiation Using MiRaGE Method. Int J Mol Sci 12(12):9265-9276
}
\examples{
\dontrun{HS_conv_id()}
}

