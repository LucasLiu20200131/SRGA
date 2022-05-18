\name{rank_vis}
\alias{rank_vis}
\title{
  Visualization of SRGA.result
}
\description{
  Draw a scatter plot with ggplot2 function geom_point to visualze SRGA.result and rank the genes based on the average rank score of signature(s).
}

\usage{
  rank_vis(SRGA.result)
}

\examples{
  data("covariate",package='SRGA')
  data("exprs",package='SRGA') # makesure the input exprs is a matrix/dataframe with gene in rows and sample in columns.
  data("Sene.marker",package='SRGA')
  set.seed(1)
  select.name = sample(rownames(exprs),100)
  example.result = signature_related_mRNA(exprs,Sene.marker,covariate,select.name,scale.flag=FALSE)
  rank_vis(example.result)
}
