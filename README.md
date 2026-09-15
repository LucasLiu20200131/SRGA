# SRGA
Signature-Related Gene Analysis

For an expression matrix, genes are denoted as $G_1, G_2, \ldots, G_n$, and samples as $S_1, S_2, \ldots, S_m$.

### Gene correlation

The Pearson correlation coefficient between genes $G_i$ and $G_j$ is calculated across the $m$ samples as:

$$
CC_{ij} =
\frac{
\sum_{z=1}^{m}
(G_{iz}-\bar{G}_i)(G_{jz}-\bar{G}_j)
}{
\sqrt{\sum_{z=1}^{m}(G_{iz}-\bar{G}_i)^2}
\sqrt{\sum_{z=1}^{m}(G_{jz}-\bar{G}_j)^2}
}
$$

where $G_{iz}$ and $G_{jz}$ represent the expression levels of genes $G_i$ and $G_j$ in sample $z$, respectively, and $\bar{G}_i$ and $\bar{G}_j$ represent their mean expression levels across samples.

### Partial correlation

To account for an independent covariate, such as tumor purity, SRGA can use partial correlation. The partial correlation coefficient between $G_i$ and $G_j$, controlling for tumor purity $P$, is calculated as:

$$
PCC_{ij} =
\frac{
CC_{ij} - CC_{iP}CC_{jP}
}{
\sqrt{1-CC_{iP}^{2}}
\sqrt{1-CC_{jP}^{2}}
}
$$

where $CC_{iP}$ and $CC_{jP}$ denote the correlation coefficients between tumor purity and the expression levels of $G_i$ and $G_j$, respectively.

Hereafter, the correlation coefficient between $G_i$ and $G_j$ is denoted as $\operatorname{cor}_{ij}$, with the corresponding statistical significance denoted as $p_{ij}$.

### Relative score

For each selected gene $G_i$, its association with another gene $G_j$ is quantified using a relative score:

$$
RS_{ij}
=
-\log_{10}(p_{ij})
\times
\operatorname{sign}(\operatorname{cor}_{ij})
$$

Genes are then ranked according to $RS_{ij}$ after removing self-correlations and infinite values.

### Signature enrichment score

For each selected gene, its ranked $RS$ gene list is subjected to GSEA against the input signatures. For signature $i$ and selected gene $G_j$, the signature value is calculated as:

$$
\operatorname{sigValue}_{ij}
=
-\log_{10}(p_{ij}^{\mathrm{GSEA}})
\times
NES_{ij}
$$

where $p_{ij}^{\mathrm{GSEA}}$ is the statistical significance of the enrichment and $NES_{ij}$ is the corresponding normalized enrichment score.

### Relative rank score

To compare enrichment results across signatures, the signature values are rescaled within each signature. For signature $i$ and gene $G_j$:

$$
RRS_{ij}
=
\frac{
\operatorname{sigValue}_{ij}
-
\min_k(\operatorname{sigValue}_{ik})
}{
\max_k(\operatorname{sigValue}_{ik})
-
\min_k(\operatorname{sigValue}_{ik})
}
$$

where the minimum and maximum are calculated across all evaluated genes $k$ for signature $i$.

The final score of gene $G_j$ is calculated as the mean relative rank score across all $n$ input signatures:

$$
\operatorname{Rank}(G_j)
=
\frac{1}{n}
\sum_{i=1}^{n} RRS_{ij}
$$


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
