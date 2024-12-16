######################################
#Functions for the update of ADMM
#######################################
##take the larger one compared to 0
pl <-function(TC)
{
  TC[TC<0]=0
  return(TC)
}

#The S/C function in paper
S_c <- function(llambda,oomega)
{
  B_a<-matrix(1,ncol=ncol(oomega),nrow=nrow(oomega))
  if(length( which(oomega > llambda) )>0 )  {B_a[oomega > llambda] <- oomega[oomega > llambda] - llambda}
  if(length( which(oomega < 0) )>0 )  {B_a[oomega < 0] <- oomega[oomega < 0]}
  if(length( which(oomega >= 0 & oomega <= llambda) )>0 )  {B_a[oomega >= 0 & oomega <= llambda] <- 0}
  return(B_a)}

#The S/b function in paper
S_b <- function(llambda,oomega)
{
  B_b<-matrix(1,ncol=ncol(oomega),nrow=nrow(oomega))
  if(length( which(oomega >= 0) )>0 ) {B_b[oomega >= 0] <- oomega[oomega >=0] - llambda}
  if(length( which(oomega < 0) )>0 ) {B_b[oomega < 0] <- oomega[oomega < 0] + llambda}
  return(B_b)
}

#The S/d function in paper
S_d <- function(llambda,oomega)
{
  B_d<-matrix(1,ncol=ncol(oomega),nrow=nrow(oomega))
  if(length( which(oomega > llambda) )>0 ) {B_d[oomega > llambda] <- oomega[oomega > llambda] - llambda}
  if(length( which(oomega < -llambda) )>0 ) { B_d[oomega < -llambda] <- oomega[oomega < -llambda] + llambda}
  if(length( which(oomega >= -llambda & oomega <= llambda) )>0 ) {B_d[oomega >= -llambda & oomega <= llambda] <- 0}
  return(B_d)
}

#The g_c function for b
G_c <- function(glambda,gomega)
{ 
  B_a<-matrix(1,ncol=ncol(gomega),nrow=nrow(gomega))
  B_a[gomega > 0] <- gomega[gomega > 0]-glambda
  B_a[gomega <= 0] <- gomega[gomega <= 0]
  return(B_a)
}

ka_z<-function(z=tz_,nnz=0.01)
{
  z[abs(z)<=nnz]=0
  return(z)
}

#The T function for b
T_c <- function(tlambda, tomega)
{
  sign(tomega)*pl(abs(tomega)-tlambda)
}

#Update of the block t for GKF-PUAL
block_update<-function(i,XX_bl,lambda,v,mu,D_star,alpha_star)
{ 
  XX_b<-as.matrix(XX_bl[, ((i-1)*(m+1)+1):(i*(m+1))  ])
  max( 1-lambda/(2*mu*norm(t(XX_b)%*%(D_star%*%alpha_star+v/mu))),0)*t(XX_b)%*%(D_star%*%alpha_star+v/mu)
}

###############################################################
#The Construction of the design matrix for the update of t
##############################################################
design_non<-function(con_s,v_a_1)
{ 
  m=length(v_a_1)-1
  mm<-0.5*m*(m+1)+m
  output_ini<-rep(0,times=mm)
  output_ini[v_a_1[con_s]]<-1
  output_ini
}

design_final<-function(XX_114514)
{
  ii<-nrow(XX_114514)
  iii<-matrix(1:ii,ncol=1)
  apply(iii,MARGIN=1,design_non_overlap,v=XX_114514)
}

XU_check<-function(ZT_,x9_=0.75)
{
  if(ZT_<x9_)
  {return(1)}
  else
    return(0)
}

design_format<-function(XX_4)
{ 
  m<-ncol(XX_4)
  n_rrow<-m+1
  n_ccol<-0.5*m*(m+1)+m
  n_re=ncol(XX_4)
  XX_re_all<-matrix(0,nrow=0,ncol=n_ccol)
  for(i_re in 1:n_re)
  { 
    XX_re<-as.vector(XX_4[,i_re])
    XX_re_1<-t(matrix(XX_re ,nrow=n_ccol,ncol=n_rrow,byrow = FALSE))
    XX_re_all<-rbind(XX_re_all,XX_re_1)
  }
  XX_re_all
}

dui_rgl<-function(a,b)
{
  resp<-sign(a)*max((abs(a)-b),0 )
  return(resp)
}

dui_rgl2<-function(a1,b)
{ 
  resp2<-c()
  resp2<-apply(a1,MARGIN=1,FUN=dui_rgl,b=b)
  resp4<-matrix(resp2,ncol=1)
  return(resp4)
}


########################################################################
##Functions to change the form of the model matching to the linear case
########################################################################
op_convert<-function(j,XX)
{ 
  m=length(XX)
  m_global=0.5*(m^2+m)
  XX_convert<-rep(0,times=m_global)
  x_1<-seq(1:(j-1))
  j_111<-rep(j,times=j-1)
  first_part<-0.5*(2*m-x_1+2)*(x_1-1)+(j_111-x_1)+rep(1,times=j-1)
  second_part<-seq(from=0.5*(2*m-j+2)*(j-1)+1, to=0.5*(2*m-j+2)*(j-1)+m-j+1)
  x_row_index<-c(first_part,second_part)
  XX_convert[x_row_index]<-XX
  XX_convert
}

fin_convert<-function(j_11,XX_2)
{ 
  XX_1<-XX_2[j_11,]
  j_1 <-matrix(seq(1:length(XX_1)),ncol=1)
  apply(j_1, MARGIN=1, FUN=op_convert,XX=XX_1)
}

fin_overall_convert<-function(XX_3)
  {
  jjjjj_over<-matrix(seq(1:nrow(XX_3)),ncol=1)
  result_conv<-apply(jjjjj_over, MARGIN=1, FUN=fin_convert,XX_2=XX_3)
  as.matrix(result_conv)
}

re_format<-function(XX_4,n_ccol)
{ n_rrow<-0.5*(n_ccol^2+n_ccol)
  n_re=ncol(XX_4)
  XX_re_all<-matrix(0,nrow=0,ncol=n_rrow+n_ccol)
  for(i_re in 1:n_re)
  { 
    XX_re<-as.vector(XX_4[,i_re])
    XX_re_1<-t(matrix(XX_re ,nrow=n_rrow,ncol=n_ccol,byrow = FALSE))
    XX_re_1<-cbind(XX_re_1,diag(n_ccol))
    XX_re_all<-rbind(XX_re_all,XX_re_1)
  }
  XX_re_all
}

re_S<-function(v_index,m,XX,x)
{ 
  XXX_1<-as.matrix(XX[(m*(v_index-1)+1):(m*v_index),])
  XXX_2<-matrix(x[v_index,],ncol=1)
  0.5*( t(XXX_1)%*%XXX_2+t(cbind(matrix(0,nrow=1,ncol=0.5*(m^2+m)),t(XXX_2) )))
}

design_non_overlap<-function(i,v)
{ 
  zet_<-nrow(v)+1
  v_i<-matrix(v[i,],nrow=1)
  v_a<-which(v_i==1) 
  index_p<-matrix(1:zet_,ncol=1)
  apply(index_p,MARGIN=1,FUN=design_non,v_a_1=v_a)
}

matrix_0<-function(XXX)
{ 
  mm<-nrow(XXX)
  nm<-ncol(XXX)
  rbind(cbind(XXX,matrix(0,nrow=mm,ncol=1)  ), matrix(0,ncol=nm+1,nrow=1)   )}

matrix_1<-function(vv)
{cbind(vv,matrix(1,nrow=nrow(vv),ncol=1))}

vector_0<-function(vv)
{
  rbind(vv,0)
}

###########################################################
#Change the form of model parameters back to the initial one
###########################################################
cons_ver<-function(vexc,fex,desig)
{ 
  fex_abs<-abs(fex)
  index_00<-which(fex_abs<0.000001)
  if(length(index_00)>0)
  {  
    x_z1<-as.matrix(desig[index_00,])
    azd1<-matrix(1:nrow(x_z1),ncol=1)
    index_11<-apply(azd1,FUN = find_111,MARGIN = 1,De_d=x_z1)
    vexc[index_11,1]=0
    vexc
   }
  else
  vexc
}

ZRTD <- function(BUP_,CON_BUP,TEX_,n_32,m_=m+1)
{   
  if(length(which(abs(TEX_)<0.001))<(n_32*m_+1))
      {c(BUP_,CON_BUP)}
}


find_111<-function(a__1,De_d)
{
  xx1<-De_d[a__1,]
  which(xx1==1)

}




