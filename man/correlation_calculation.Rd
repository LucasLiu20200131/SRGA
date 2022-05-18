\name{correlation_calculation}
\alias{correlation_calculation}
\title{
  function to do correlation analysis
}
\description{
  function to do correlation analysis, if covariate=NULL (default), use rcorr function from Hmisc package, otherwise partial correlation analysis performed. for more information, please refer to signature_related_mRNA function
}
\usage{
  correlation_calculation(mRNA_exp,covariate,select.name)
}


\examples{
  # cor
  data("exprs",package='SRGA')
  select.name = sample(rownames(exprs),100)
  RS=correlation_calculation(mRNA_exp=exprs,select.name=select.name)
  
  # PCC
  data("covariate",package='SRGA')
  data("exprs",package='SRGA')
  set.seed(1)
  select.name = sample(rownames(exprs),100)
  RS=correlation_calculation(exprs,covariate,select.name)
}
