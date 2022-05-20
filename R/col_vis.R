col_vis <- function(SRGA.result,log2.flag=TRUE){
  for(pkg in c('tidyverse')){
    if(!requireNamespace(pkg,quietly=T)){
      stop(paste("The ",pkg," package needed for this function to work. Please install it.",sep=""),
           call. = FALSE)
    }else{
      library(pkg,warn.conflicts=F,character.only = T)
    }
  }
  all_pathway=SRGA.result %>%
    filter(pval <0.05)%>%
    group_by(pathway)%>%
    count()%>%
    arrange(n)
  if(log2.flag){
    n_pathway = sum(all_pathway$n>1)
    message(paste0('There are ',n_pathway,' Signature(s) with number of related-genes more than 0 after log2 transforming, plotting them..'))
    all_pathway = all_pathway[all_pathway$n!=1,]
    all_pathway$pathway = factor(all_pathway$pathway,levels= all_pathway$pathway)
    p1=ggplot(all_pathway)+
      geom_col(aes(x=pathway,y=log2(n)),fill='#61c8c0',color='black')+
      coord_flip()+
      theme_classic()+
      theme(axis.ticks.y = element_blank())+
      labs(x=NULL,y='log2(Number of signature-related genes)')
  }else{
    n_pathway = sum(all_pathway$n>0)
    message(paste0('There are ',n_pathway,' Signature(s) with number of related-genes more than 0, plotting them..'))
    all_pathway$pathway = factor(all_pathway$pathway,levels= all_pathway$pathway)
    p1=ggplot(all_pathway)+
      geom_col(aes(x=pathway,y=n),fill='#61c8c0',color='black')+
      coord_flip()+
      theme_classic()+
      theme(axis.ticks.y = element_blank())+
      labs(x=NULL,y='Number of signature-related genes')
  }
  return(p1)
}
