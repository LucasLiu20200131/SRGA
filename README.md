# SRGA
signature-related gene analysis

# Install
```devtools::install_github("LucasLiu20200131/SRGA")```

# Usage
```
data("covariate",package='SRGA')
data("exprs",package='SRGA') # makesure the input exprs is a matrix/dataframe with gene in rows and sample in columns.
data("Sene.marker",package='SRGA')
select.name = sample(rownames(exprs),100)
example.result = signature_related_mRNA(exprs,Sene.marker,covariate,select.name,scale.flag=FALSE)
}
```
