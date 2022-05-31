\name{gene_exclude}
\alias{gene_exclude}
\title{
  Pre-process of gene expression matrix
}
\description{
  Remove genes (rows) with zero expression in more than 'ex.per' samples (columns).
}

\usage{
gene_exclude(exprs,ex.per=0.3)
}

\arguments{
  \item{exprs}{Gene expression matrix}
  \item{ex.per}{Gene with this proportion of zero expression in samples that to be excluded}
}

\examples{
data("exprs",package='SRGA')
exprs = gene_exclude(exprs,ex.per=0.3)
}
