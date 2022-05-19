net_vis <- function(SRGA.result,top.gene=5){
  for(pkg in c('igraph')){
    if(!requireNamespace(pkg,quietly=T)){
      stop(paste("The ",pkg," package needed for this function to work. Please install it.",sep=""),
           call. = FALSE)
    }else{
      library(pkg,warn.conflicts=F,character.only = T)
    }
  }
  top_df = data.frame()
  for (p in unique(SRGA.result$pathway)){
    p.se = SRGA.result[SRGA.result$pathway==p,]
    p.se <- p.se[order(p.se$sigValue,decreasing = T),]
    se = p.se[1:top.gene,c('pathway','mRNA')]
    top_df = rbind(top_df,se)
  }
  g1 = graph_from_data_frame(top_df,directed = F)
  V(g1)$color[V(g1)$name %in% unique(top_df$pathway)]='#fd3737'
  V(g1)$color[!V(g1)$name %in% unique(top_df$pathway)]='#a3bae2'
  V(g1)$size[V(g1)$name %in% unique(top_df$pathway)]=8
  V(g1)$size[!V(g1)$name %in% unique(top_df$pathway)]=4
  plot(g1,
       layout=layout.auto,
       vertex.label=V(g1)$name,
       vertex.label.cex=ifelse(V(g1)$name %in% unique(top_df$pathway),0.8,0.4),
       vertex.label.dist=ifelse(V(g1)$name %in% unique(top_df$pathway),0,0.9),
       vertex.label.color='black')
  legend('topright',legend=c('Signature','Gene'),pt.cex= 1.5,
         pch=16,col= c('#fd3737','#a3bae2'),pt.bg = 'white',title.col = 'white',
         border = 'white',x.intersp = 0.5,text.col='black',cex = 0.8)
  return(top_df)
}
