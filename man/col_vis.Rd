\name{col_vis}
\alias{col_vis}
\title{
  Visualization of SRGA.result
}
\description{
  Draw a bar plot with ggplot2 function geom_col to visualze SRGA.result and count the number of related genes for each signature.
}

\usage{
  col_vis(SRGA.result,log2.flag=TRUE)
}

\arguments{
  \item{SRGA.result}{A data.frame object, output of signature_related_mRNA()}
  \item{log2.flag}{Should the number of signature-related genes transformed by log2 when visualiztion? Defualt TRUE}
}

\examples{
  data("covariate",package='SRGA')
  data("exprs",package='SRGA') # makesure the input exprs is a matrix/dataframe with gene in rows and sample in columns.
  exprs = gene_exclude(exprs,ex.per=0.3)
  data("Sene.marker",package='SRGA')
  set.seed(1)
  select.name = sample(rownames(exprs),100)
  example.result = signature_related_mRNA(exprs,Sene.marker,covariate,select.name,scale.flag=FALSE)
  col_vis(example.result,log2.flag = T)
  col_vis(example.result,log2.flag = F)
}
