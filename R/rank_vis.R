rank_vis <- function(SRGA.result,top.show=5){
  for(pkg in c('scales','ggplot2')){
    if(!requireNamespace(pkg,quietly=T)){
      stop(paste("The ",pkg," package needed for this function to work. Please install it.",sep=""),
           call. = FALSE)
    }else{
      library(pkg,warn.conflicts=F,character.only = T)
    }
  }
  pathway_rank_score=lapply(unique(SRGA.result$pathway),function(x)
  { p = SRGA.result[SRGA.result$pathway==x,]
    if(sum(is.na(p$sigValue))==nrow(p)){
      score = rep(0,nrow(p))
    }else{
    r = p[order(p$sigValue,na.last = F),]
    r$score = scales::rescale(1:nrow(r),c(0,1))
    score  = r$score[match(p$mRNA,r$mRNA)]}
    names(score) = p$mRNA
    return(score)
  })
  names(pathway_rank_score) = unique(SRGA.result$pathway)
  all_rank_score <- Reduce(function(x,y) cbind(x,y),pathway_rank_score)
  colnames(all_rank_score) = names(pathway_rank_score)
  all_rank_score = data.frame(all_rank_score,check.names = F)
  all_rank_score$rank <- rowMeans(all_rank_score)
  all_rank_score <- all_rank_score[order(all_rank_score$rank,decreasing = T),]
  all_rank_score$number = seq(1:nrow(all_rank_score))
  p1=ggplot(data=all_rank_score)+
    geom_point(aes(x=number,y=rank),color='#61c8c0')+
    theme_classic()+
    xlab('Number of genes')+
    ylab('Rank')
  if(top.show){
    gene.show = paste(rownames(all_rank_score)[1:top.show],collapse = '\n')
    p1=p1+
      geom_segment(aes(x = number[1], y = rank[1], xend = nrow(all_rank_score)/10, yend = 0.3),colour = "black",lty = 1)+
      geom_label(x=nrow(all_rank_score)/10,y=0.3,label=gene.show,label.size=0.5)
  }
  return(p1)
}
