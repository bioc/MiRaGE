MiRaGE <- function(gene_exp,location="local",species="MM",ID="refseq",method="mean",test="ks",conv="conserved",species_force=T,ID_force=T,conv_force=T)
{
tb1 <- getMiRaGEData(location=location, species=species, ID=ID,method=method,test=test,conv=conv,species_force=species_force,ID_force=ID_force,conv_force=conv_force)
return(MiRaGEAnalysis(gene_exp, tb1, method=method, test=test))
}
