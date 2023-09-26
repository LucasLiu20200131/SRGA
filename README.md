# SRGA
Signature-Related Gene Analysis

Briefly, expression matrix was obtained with gene in row and sample in column. After quality control step was applied to remove unnecessary genes with high percentage of zero expression in samples that might confound results, the correlations among genes were calculated and genes were ranked based on a scoring formula. Then, the ranked gene list and interested signatures were subjected to GSEA method to get a score value for each gene within each signature. Finally, the relative rank of genes was determined based on average signature scores.
In the correlation calculation step, Pearson correlation analysis was performed. For an expression matrix, genes defined as G1, G2, …, Gn, and samples defined as S1, S2, …, Sm. The correlation coefficient CCij between Gi and Gj (i, j  n) was calculated as follow:
〖CC〗_ij " = " z=1n(Giz- Gi)(Gjz- Gj)∑_"z=1" ^"n" ▒("G" _"iz"  "-" (G_"i"  ) ̅ )^"2"  z=1nGjz-Gj2
Considering the effect of other independent factor, like in this study, tumor purity, we introduce partial correlation analysis to alleviate the impact of covariate. The partial correlation coefficient PCCij was calculated as follow:
〖PCC〗_ij " = "  (〖CC〗_"ij" -〖CC〗_"ip" *〖CC〗_"jp" )⁄(√(1-〖〖CC〗_"ip" 〗^2 )*√(1-〖〖CC〗_"jp" 〗^2 ))
where CCip and CCjp indicated the correlation coefficients between the expression of Gi and tumor purity, and the expression of Gj and tumor purity, respectively. The correlation coefficient between Gi and Gj was represented as cor(ij) hereafter, and we also obtained the p value and defined it as p(ij). For each gene from the cluster markers, we calculated its relative score (RS) with the rest of other genes as follow:
RS(ij)"= "-log10(p(ij))*sign(cor(ij))
	In the GSEA-based signature modulators analysis step, we ranked all genes based on their relative scores after removing self-correlation and infinite values. For a given gene, the RS rank list and 13 senescent signatures, including 12 aging-related pathways from AingAtlas and the Sene.marker gene set defined previously in this study, were feed into fgsea function. We calculated a sigValue for every selected gene in each signature as follow:
sigValue "= "-log10(p(i))*NES(i)
	Where p was the statistical significance of the enrichment score and NES was the normalized enrichment score for all dataset permutations.
	Unlike ImmLnc, we considered the relative position in signatures instead of absolute value to select significant genes. For a given signature i (i n, n = the total number of signatures input), we assigned each gene Gj a relative rank score (RRS) as follow:
RRS(ij)"= "  (sigValue(G_j)-min⁡(i)  )⁄( max⁡(i) )-min⁡(i)
	And the final rank of genes was calculated as follow:
Rank(Gj)"= "  (∑_(i=1)^n▒〖RRS(ij)〗  )⁄( n)

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
