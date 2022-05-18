rank_vis <- function(SRGA.result){
  for(pkg in c('scales','ggplot2')){
    if(!requireNamespace(pkg,quietly=T)){
      stop(paste("The ",pkg," package needed for this function to work. Please install it.",sep=""),
           call. = FALSE)
    }else{
      library(pkg,warn.conflicts=F,character.only = T)
    }
  }
  pathway_rank_score=lapply(unique(SRGA.result$pathway),function(x)
  {
    p = SRGA.result[SRGA.result$pathway==x,]
    r = p[order(p$sigValue),]
    r$score = scales::rescale(1:nrow(r),c(0,1))
    score  = r$score[match(p$mRNA,r$mRNA)]
    names(score) = p$mRNA
    return(score)
  })
  names(pathway_rank_score) = unique(SRGA.result$pathway)
  all_rank_score <- Reduce(function(x,y) cbind(x,y),pathway_rank_score)
  colnames(all_rank_score) = names(pathway_rank_score)
  all_rank_score = data.frame(all_rank_score,check.names = F)
  all_rank_score$rank <- rowMeans(all_rank_score)
  all_rank_score <- all_rank_score[order(all_rank_score$rank,decreasing = T),]
  top20 <- rownames(all_rank_score)[1:20]
  all_rank_score$number = seq(1:nrow(all_rank_score))
  p1=ggplot(data=all_rank_score)+
    geom_point(aes(x=number,y=rank),color='#61c8c0')+
    theme_classic()+
    xlab('Number of genes')+
    ylab('Rank')
  return(p1)
}