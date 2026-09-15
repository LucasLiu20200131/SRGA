# SRGA

**Signature Related Gene Analysis**

SRGA is an R package for identifying genes associated with predefined biological signatures based on gene expression correlation and gene set enrichment analysis.

For an expression matrix, genes are denoted as $G_1, G_2, \ldots, G_n$, and samples are denoted as $S_1, S_2, \ldots, S_m$.

## Method

### Gene correlation

The Pearson correlation coefficient between genes $G_i$ and $G_j$ is calculated across the $m$ samples as:

```math
CC_{ij} =
\frac{
\sum_{z=1}^{m}
(G_{iz}-\bar{G}_i)(G_{jz}-\bar{G}_j)
}{
\sqrt{\sum_{z=1}^{m}(G_{iz}-\bar{G}_i)^2}
\sqrt{\sum_{z=1}^{m}(G_{jz}-\bar{G}_j)^2}
}
```

where $G_{iz}$ and $G_{jz}$ represent the expression levels of genes $G_i$ and $G_j$ in sample $z$, respectively, and $\bar{G}_i$ and $\bar{G}_j$ represent their mean expression levels across samples.

### Partial correlation

To account for a potential covariate, such as tumor purity, SRGA can alternatively use partial correlation. For a covariate $P$, the partial correlation coefficient between $G_i$ and $G_j$ is calculated as:

```math
PCC_{ij} =
\frac{
CC_{ij} - CC_{iP}CC_{jP}
}{
\sqrt{1-CC_{iP}^{2}}
\sqrt{1-CC_{jP}^{2}}
}
```

where $CC_{iP}$ and $CC_{jP}$ denote the correlation coefficients between the covariate $P$ and the expression levels of $G_i$ and $G_j$, respectively.

Hereafter, $r_{ij}$ denotes the correlation coefficient used for genes $G_i$ and $G_j$, and $p_{ij}$ denotes its corresponding P value.

### Relative score

For each selected gene $G_i$, its association with another gene $G_j$ is quantified using a relative score, $RS_{ij}$:

```math
RS_{ij} =
-\log_{10}(p_{ij})
\times
\mathrm{sign}(r_{ij})
```

This formulation incorporates both the statistical evidence and direction of the correlation. Genes are ranked according to their relative scores after self correlations and infinite values are removed.

### Signature enrichment score

For each selected gene, the ranked gene list derived from the relative scores is subjected to gene set enrichment analysis against the input signatures.

For signature $i$ and selected gene $G_j$, the signature enrichment value $SV_{ij}$ is calculated as:

```math
SV_{ij} =
-\log_{10}(p_{ij}^{\mathrm{GSEA}})
\times
NES_{ij}
```

where $p_{ij}^{\mathrm{GSEA}}$ denotes the GSEA P value and $NES_{ij}$ denotes the corresponding normalized enrichment score.

### Relative rank score

To facilitate comparison across signatures, the signature enrichment values are rescaled within each signature.

For signature $i$ and selected gene $G_j$, the relative rank score $RRS_{ij}$ is calculated as:

```math
RRS_{ij} =
\frac{
SV_{ij} - \min_k(SV_{ik})
}{
\max_k(SV_{ik}) - \min_k(SV_{ik})
}
```

where the minimum and maximum are calculated across all evaluated genes $k$ for signature $i$.

If $K$ signatures are provided, the final score of gene $G_j$ is calculated as the mean relative rank score across all signatures:

```math
\mathrm{Rank}(G_j) =
\frac{1}{K}
\sum_{i=1}^{K} RRS_{ij}
```

Genes can subsequently be prioritized according to this final score, with higher scores indicating stronger overall associations with the input signatures.

# Installation

Install the required packages and SRGA from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("fgsea")

install.packages(
    c("scales", "Hmisc", "igraph", "tidyverse", "devtools")
)

devtools::install_github("LucasLiu20200131/SRGA")
```

# Usage

## 1. Load SRGA and example data

```r
library(SRGA)

data("covariate", package = "SRGA")
data("exprs", package = "SRGA")
data("Sene.marker", package = "SRGA")
```

The input expression data should be provided as a matrix or data frame with genes in rows and samples in columns.

Genes with a high proportion of zero expression values can be removed using `gene_exclude()`:

```r
exprs <- gene_exclude(
    exprs,
    ex.per = 0.3
)
```

For demonstration, randomly select 100 genes:

```r
set.seed(1)

select.name <- sample(
    rownames(exprs),
    100
)
```

## 2. Run SRGA

```r
example.result <- signature_related_mRNA(
    exprs,
    Sene.marker,
    covariate,
    select.name,
    scale.flag = FALSE
)
```

## 3. Visualize the results

### 3.1 Signature associated gene counts

Use `col_vis()` to visualize the number of associated genes for each signature.

With log2 transformation:

```r
col_vis(
    example.result,
    log2.flag = TRUE
)
```

![col\_vis\_1](https://github.com/LucasLiu20200131/images/blob/main/git_image/col_vis_1.png)

Without log2 transformation:

```r
col_vis(
    example.result,
    log2.flag = FALSE
)
```

![col\_vis\_2](https://github.com/LucasLiu20200131/images/blob/main/git_image/col_vis_2.png)

### 3.2 Gene ranking

Use `rank_vis()` to visualize genes according to their average relative rank scores across signatures:

```r
rank_vis(example.result)
```

![rank\_vis\_1](https://github.com/LucasLiu20200131/images/blob/main/git_image/rank_vis_1.png)

Detailed ranking information can also be returned:

```r
rank_info <- rank_vis(
    example.result,
    res.return = "rank"
)
```

### 3.3 Signature gene network

Use `net_vis()` to visualize the top genes associated with each signature. The function also returns information on the corresponding signature gene pairs:

```r
net_info <- net_vis(
    example.result,
    top.gene = 10
)
```

![net\_vis\_1](https://github.com/LucasLiu20200131/images/blob/main/git_image/net_vis_1.png)
