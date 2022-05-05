signature_related_mRNA <- function(mRNA_exp,signatures,covariate,select.name,scale.flag=FALSE){
  RS=correlation_calculation(mRNA_exp,covariate,select.name)
  fgsea_all=fgsea_calculation(RS,signatures,scale.flag)
  return(fgsea_all)
}
