gene_exclude <- function(exprs,ex.per = 0.3){
  mRow_exclude <- apply(exprs,1,function(x) {return(sum(x==0)/length((x))>=ex.per)})
  return(exprs[!mRow_exclude,])
}
