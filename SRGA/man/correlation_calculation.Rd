\name{correlation_calculation}
\alias{correlation_calculation}
\title{
function to do partial correlation analysis
}
\description{
please refer to signature_related_mRNA function
}
\usage{
correlation_calculation(mRNA_exp,covariate,select.name)
}


\examples{
data("covariate",package='SRGA')
data("exprs",package='SRGA')
select.name = sample(rownames(exprs),100)
RS=correlation_calculation(exprs,covariate,select.name)
