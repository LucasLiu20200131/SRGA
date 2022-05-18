\name{net_vis}
\alias{net_vis}
\title{
  Visualization of SRGA.result
}
\description{
  Draw a network with package igraph to visualze SRGA.result and display the top genes related to each signature.
}

\usage{
  net_vis(SRGA.result,top.gene=5)
}

\arguments{
  \item{SRGA.result}{A data.frame object, output of signature_related_mRNA()}
  \item{top.gene}{Number of top genes related to each signature that to be displayed. Default 5.}
}

\examples{
  data("covariate",package='SRGA')
  data("exprs",package='SRGA') # makesure the input exprs is a matrix/dataframe with gene in rows and sample in columns.
  data("Sene.marker",package='SRGA')
  set.seed(1)
  select.name = sample(rownames(exprs),100)
  example.result = signature_related_mRNA(exprs,Sene.marker,covariate,select.name,scale.flag=FALSE)
  net_info = net_vis(example.result,top.gene = 10)
}
