\name{signature_related_mRNA}
\alias{signature_related_mRNA}
\title{
  function to identify signature-related genes
}
\description{
  A warpper first perform correlation analysis for each gene based on expressions, then do GSEA analysis to find out signature-related genes
  We recommond set scale.flag=TRUE if only one signature input, otherwise you should rank these genes based on their sigValue in each signature and take average rank into consideration so that the most relavant gene could be selected.
  We also provide three visualization methods, please refer to examples and individual function page for details.
}

\usage{
  signature_related_mRNA(mRNA_exp, signatures, covariate, select.name, scale.flag = FALSE)
}

\arguments{
  \item{mRNA_exp}{a data.frame object with genes in row and samples in column.}
  \item{signatures}{a list object input for fgsea. Could be pathways from GO, KEGG, Reactome, or gene sets from MSigDB database, or gene families, or manual curated signatures}
  \item{covariate}{A vector of covariate needed to be excluded when calculating the linear relationship between two genes. If NULL, perform pearson correlation analysis using Hmisc::rcorr(), otherwise perform partial correlation analysis. Default NULL.}
  \item{select.name}{A vector contains the interested gene names. If no preferred genes, set this parameter as 'rownames(mRNA_exp)'.}
  \item{scale.flag}{Whether rescale the Sigvalue of each gene to range between -1 and 1? Default FALSE.}
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
  rank_vis(example.result)
  net_info = net_vis(example.result,top.gene = 10)
}
