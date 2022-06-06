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
First, attach this package and import data.
```
library("SRGA")
data("covariate",package='SRGA')
data("exprs",package='SRGA') # makesure the input exprs is a matrix/dataframe with gene in rows and sample in columns.
exprs = gene_exclude(exprs,ex.per=0.3)
data("Sene.marker",package='SRGA')
set.seed(1)
select.name = sample(rownames(exprs),100)
```
Second, run SRGA
```
example.result = signature_related_mRNA(exprs,Sene.marker,covariate,select.name,scale.flag=FALSE)
```
We can now visualize results with 3 ways.  
1.1 Draw a bar plot and count the number of related genes for each signature. Below are log2 normalized.
```
col_vis(example.result,log2.flag = T)
```
![col_vis_1](https://github.com/LucasLiu20200131/images/blob/main/git_image/col_vis_1.png)  
1.2 number are not log2 normalized.
```
col_vis(example.result,log2.flag = F)
```
![col_vis_2](https://github.com/LucasLiu20200131/images/blob/main/git_image/col_vis_2.png)  
2.1 Draw a scatter plot and rank the genes based on the average rank score of signature(s).
```
rank_vis(example.result)
```
![rank_vis_1](https://github.com/LucasLiu20200131/images/blob/main/git_image/rank_vis_1.png)  
we can also obtained the detailed rank information with:
```
rank_info = rank_vis(example.result,res.return="rank")
```
3.1 Draw a network and display the top genes related to each signature. Also return the signature-gene pairs information.
```
net_info = net_vis(example.result,top.gene = 10)
```
![net_vis_1](https://github.com/LucasLiu20200131/images/blob/main/git_image/net_vis_1.png)
