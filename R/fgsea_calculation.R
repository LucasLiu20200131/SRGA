fgsea_calculation <-function (RS,signatures,scale.flag=FALSE){
  for(pkg in c('fgsea','scales')){
    if(!requireNamespace(pkg,quietly=T)){
      stop(paste("The ",pkg," package needed for this function to work. Please install it.",sep=""),
           call. = FALSE)
    }else{
      library(pkg,warn.conflicts=F,character.only = T)
    }
  }
  fgsea_all <- c()
  for(i in 1:nrow(RS)){
    ##remove self-correlation
    mRNA <- rownames(RS)[i]
    ranks <- RS[i,-grep(mRNA,colnames(RS))]
    ##remove Infs
    ranks <- ranks[!is.infinite(ranks)]
    ranks <- sapply(ranks, as.numeric)
    fgseaRes <- fgsea(signatures, ranks, minSize=1, maxSize=5000, nperm=10000)
    addname = names(signatures)[!names(signatures) %in% fgseaRes$pathway]
    if(length(addname)){
    addrow = data.frame(pathway=addname,
                        pval=NA,
                        padj=NA,
                        ES=NA,
                        NES=NA,
                        nMoreExtreme=NA,
                        size=NA,
                        leadingEdge=NA)
    fgseaRes = rbind(fgseaRes,addrow)}
    fgseaRes$sigValue <- -log10(fgseaRes$padj) * fgseaRes$NES
    fgseaRes <- cbind(mRNA,fgseaRes)
    fgsea_all <- rbind(fgsea_all,fgseaRes)
  }
  if(scale.flag){
    fgsea_all$sigValue <- scales::rescale(fgsea_all$sigValue,c(-1,1))
  }
  return(fgsea_all)
}
