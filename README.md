# SRGA
Signature-Related Gene Analysis

# Install
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("fgsea")
install.packages(c("scales","Hmisc","igraph","tidyverse")) 
devtools::install_github("LucasLiu20200131/SRGA")
```

# Usage
```
library("SRGA")
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
```
