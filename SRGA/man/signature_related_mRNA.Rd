\name{signature_related_mRNA}
\alias{signature_related_mRNA}
\title{
function to identify signature-related genes
}
\description{
    A warpper to do first partial correlation analysis among expressions, then do GSEA analysis to find out signature-related genes
    We recommond set scale.flag=TRUE if only one signature input, otherwise you should rank these genes based on their sigValue in each signature and take average rank into consideration so that the most relavant gene could be selected.
}
\usage{
signature_related_mRNA(mRNA_exp, signatures, covariate, select.name, scale.flag = FALSE)
}

\examples{
data("covariate",package='SRGA')
data("exprs",package='SRGA')
data("Sene.marker",package='SRGA')
select.name = sample(rownames(exprs),100)
example.result = signature_related_mRNA(exprs,Sene.marker,covariate,select.name,scale.flag=FALSE)
}
