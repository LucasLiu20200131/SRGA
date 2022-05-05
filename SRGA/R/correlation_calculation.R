correlation_calculation <- function(mRNA_exp,covariate,select.name){
  fun_mtx_pcr <- function(x,y,z){
    r12=cor(t(x),t(y))
    r13=cor(t(x),z)
    r23=cor(z,t(y))
    r123=r13%*%r23
    rup=r12-r123
    rd1=sqrt(1-r13*r13)
    rd2=sqrt(1-r23*r23)
    rd=rd1%*%rd2
    rrr=rup/rd
    return(rrr)
  }
  n=length(mRNA_exp)
  gn=1
  pcor <- fun_mtx_pcr(mRNA_exp,mRNA_exp,covariate)
  statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  p.value <- 2*pnorm(-abs(statistic))
  RS <- -log10(p.value)*sign(pcor)
  rownames(RS) <- rownames(mRNA_exp)
  colnames(RS) <- rownames(mRNA_exp)
  return(RS = RS[select.name,])
}
