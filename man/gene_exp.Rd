\name{gene_exp}
\alias{gene_exp}
\docType{data}
\title{
This ExpressionSet represents gene expression of lung cancer cell line one day after let-7a transfection.
}
\description{
Gene expression of lung cancer cell line one day after let-7a transfection and negative control.
Each of them consists of two biological replicates. 
}
\usage{
data(gene_exp)
}
\format{
gene_id in featureData is mainly, RefSeq mRNA, which is supplemented by GenBank ID if RefSeq mRNA is not available. The expressions are two biological replicates: one day after transfection of let-7a and two biological replicates one day after transfection of mock (negative control), respectively. Sample_name in phenoData represents these characteristic features of samples.
}
\source{
This comes from experiments described in Takahashi  et al. (2009).
} 
\references{
Takahashi, Y., Forrest, A.A.R., Maeno, E., Hashimoto, T., Daub, C.O., Yasuda,
J.: MiR-107 and MiR-185 Can Induce Cell cycle Arrest in Human Non Small Cell
Lung Cancer Cell Lines? PLoS One 4, e6677 (2009)
}
