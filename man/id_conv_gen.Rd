\name{id_conv_gen}
\alias{id_conv_gen}
\title{
Conversion table between refseq id vs various gene id/probe id generation function
}
\description{
This function generate id_conv, which describes conversion between refseq gene id vs various gene id/prob id, from the scratch}
\usage{
id_conv_gen(SP="MM",gene_id="all",destdir="./")
}
\arguments{
  \item{SP}{Species for which gene id conversion table is generated. "MM" for mouse, "HS" for human. Default is "MM"}
  \item{gene_id}{gene id for which conversion table toward refseq gene is generated. If "all" is specified, tables for all of available gene ids are generated. Default is "all".} 
  \item{destdir}{Directory where data sets are saved. Default is "./".}
 }
\details{
This function generates gene id conversion table, id_conv, 
  used in the MiRaGE server described in Yoshizawa et al (2011).
 It saves data frames  in the destination directory.
The name of data set including each id_conv is the form, "SP_refseq_to_gene_id_id", where SP is either MM or HS and gene_id is gene id defined in biomaRt attributes.}
\note{
If the number of empty rows of each id_conv is less than 20000, its usage is strongly discouraged because the resulting P-values are
not trusted. 
}
\references{
Yoshizawa, M., Taguchi, Y-h., and Yasuda, J. (2011), Inference of Gene Regulation via miRNAs During ES Cell
Differentiation Using MiRaGE Method. Int J Mol Sci 12(12):9265-9276
}
\examples{
\dontrun{id_conv_gen(gene_id="ensembl_transcript_id")}
}

