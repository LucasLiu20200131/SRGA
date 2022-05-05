\name{fgsea_calculation}
\alias{fgsea_calculation}
\title{
function to do GSEA with obtained rank matrics
}
\description{
please refer to signature_related_mRNA function
}
\usage{
fgsea_calculation(RS, pathways, scale.flag)
}

\examples{
data("covariate",package='SRGA')
data("exprs",package='SRGA')
select.name = sample(rownames(exprs),100)
RS=correlation_calculation(exprs,covariate,select.name)
data("Sene.marker",package='SRGA')
fgsea_all=fgsea_calculation(RS,Sene.marker,scale.flag=FALSE)
