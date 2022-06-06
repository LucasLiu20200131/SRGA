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

\arguments{
  \item{SRGA.result}{A data.frame object, output of signature_related_mRNA()}
  \item{top.show}{Number of top genes to be displayed. Set 0 or FALSE to show nothing. Default 5.}
  \item{res.return}{return a graph in ggplot2 form or the rank information. Default 'graph', other arguments will return the rank information.}
}

\examples{
  data("covariate",package='SRGA')
  data("exprs",package='SRGA') # makesure the input exprs is a matrix/dataframe with gene in rows and sample in columns.
  data("Sene.marker",package='SRGA')
  set.seed(1)
  select.name = sample(rownames(exprs),100)
  example.result = signature_related_mRNA(exprs,Sene.marker,covariate,select.name,scale.flag=FALSE)
  rank_vis(example.result)
  rank_info = rank_vis(example.result,res.return="rank")
}
