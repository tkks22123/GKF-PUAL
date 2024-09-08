#####################
# environment
#####################
#packages
library(mvtnorm)
library(readmnist)
library(kernlab)
library(expm)
library(scatterplot3d)
library(kernlab)
library(tsne)
library(pointr)
#times
iter_times=1300
#tool functions
source('F:/ucl/pr/result/script/EX_fun_PU.R')

data111_<-read.csv("F:/ucl/pr/dataset/112/ngood/Online Retail.csv",header = T)[c(188341:188641,512708:513008),]
ran_start<-44324
set.seed(ran_start)
mg_meaningless<-rnorm(nrow(data111_), mean =0, sd = 1)
for(mgm in 1:30)
{ mg_m1<-rnorm(nrow(data111_), mean =0, sd = 1)
mg_meaningless<-cbind(mg_meaningless,mg_m1)
}
data111_<-cbind(data111_,mg_meaningless)
PP_cro<-c() 
Zero_0<-c()
error_rate_overall <- c()


ran_seed<-ran_start

for(ran_seed in (ran_start-28))
{set.seed(ran_seed)
  ranran<-sample(1:nrow(data111_))
  data1<-data111_[ranran,]
  data11<-data1[,c(1:4,6:13)]
  data12<-data1[,5]
  a2=c(50,40,-75,-10)
  data1<-data1[data12=="Germany"|data12=="United Kingdom",]
  data1<- data1[sample(1:nrow(data1)),]
  data11<-data1[,c(1:4,6:13)]
  data12<-data1[,5]
  data12[data12=="Germany"] <- -1
  data12[data12=="United Kingdom"] <- 1
  data12<-as.numeric(data12)
  defpar<-c(119,6,2.56,0.2814)
  fre_label=0.9#0.9, 2.2
  p_trte=0.7
  re_sig=0.001
  p.pi<-length(data12[data12==1])/length(data12)
  
  
  data13 <- data12
  ta<-sample(1:length(data13[data13==1] ))[1:(fre_label*length(data13[data13==1]))]
  data12[data13==1][ta]<--1

  
  data11<-data.matrix(data11, rownames.force = NA)
  data11<-data11
  
  data112<-data11
  
  data_non<-data11[data12==-1,]
  #data_non<-data_non[1:1000,]
  data_fra<-data11[data12==1,]
  data_fra<-data_fra
  
  
  Testset1<-data13[1:(p_trte*length(data13))]
  TY<-matrix(Testset1,ncol=1)
  
  # number of features
  n_f<-ncol(data11)
  Testset2<-data11[1:(p_trte*length(data13)),]
  TX<-matrix(Testset2,ncol=n_f)
  ##$$
  Testset3<-c(data13[(p_trte*length(data13)+1):length(data13)]
  )
  TY2<-matrix(Testset3,ncol=1)
  
  
  
  Testset4<-rbind(data11[(p_trte*length(data13)+1):length(data13),]
  )
  TX2<-matrix(Testset4,ncol=n_f)

###########################
# Load PU Data
###########################
# sample size, positive and negative
s_size <- length(which(TY==1))
s_size_n <-length(which(TY==-1))
# size of labeled positives
s_lp <-length(which(data12[1:(0.7*length(data13))]==1))
#No. of labeled in positives
lp <- which(data12[1:(0.7*length(data13))]==1)
lu <-intersect( which((TY==1)), which(data12[1:(0.7*length(data13))]==-1))


#arti-positive
X_P <- TX[TY==1,]
#arti-negative
X_N <- TX[TY=-1,]

##@@


#input
X_l <-  TX[data12[1:(0.7*length(data13))]==1,]
X_u <- TX[data12[1:(0.7*length(data13))]==-1,]
Y_uun1 <- matrix(1,ncol=1,nrow=s_size-s_lp)
Y_uun2 <- matrix(-1,ncol=1,nrow=s_size_n)
y_uun<- rbind(Y_uun1,Y_uun2)



#size of positive and negative
n_l <- nrow(X_l)
n_u <- nrow(X_u)
#overall sample size
nn <- n_l+n_u
X_lu<-rbind(X_l,X_u)

###### find neighbors
theta_K<-1
W<-matrix(0,nrow=nn,ncol=nn)
Phi_lu<-matrix(0,nrow=nn,ncol=nn)
sk<-sdcc(defpar)

n_32<-2
###the similar matrix W

for(j in 1:nn)
{ for(i in 1:nn)
{
  W[i,j]<-exp(- norm( matrix(X_lu[i,]-X_lu[j,],ncol = 1) ,type="2")^2 /2/theta_K)
}
 
}
diag(W)<-rep(0,times=nn)


D<-diag(rowSums(W))
L_org <- D-W
X_l_org <- X_l
X_u_org <- X_u
X_P_org <- X_P
X_N_org <- X_N
X_lu_org <- X_lu
che<-inte()
lp_org <- lp
lu_org <- lu
n_l_org <- n_l
n_u_org <- n_u
nn_org <- nn
s_size_n_org <-s_size_n
s_size_org <- s_size
s_lp_org <- s_lp
y_uun_org <- y_uun
ptr("tr3", "sk")
times_cross <- 1
as.3<-0.02
dup1<-dup_AD(TY2)
#nn <- nn_org%/%times_cross
s_size <-s_size_org 
s_size_n <-s_size_n_org 
s_lp <- s_lp
n_l <- n_l_org 
n_u <- n_u_org
nn <- n_l+n_u

#labels
y_l <- matrix(rep(1,times=n_l),ncol=1)
y_u <- matrix(rep(-1,times=n_u),ncol=1)
y_lu <- rbind(y_l,y_u)

Y_l <- matrix(0,ncol=n_l,nrow=n_l)
Y_u <- matrix(0,ncol=n_u,nrow=n_u)
Y_lu <- matrix(0,ncol=nn,nrow=nn)
lunn<-che
diag(Y_l) <- y_l
diag(Y_u) <- y_u
diag(Y_lu) <- y_lu

F_pre <- -10000000
timessss<-1
bv<-0.85
mu_2<-0.01
mu_3<-0.01
mu_4<-0.01
mu_3_star<-0.001
Phi_test<-TX2
for(lambda_1 in c(15))
for(lambda_2 in c(1))
for(C_1 in c(1)) 
for(C_2 in c(0.0264))

for(mu_1 in 0.01)
  {{{{{

  tact<-dup1
for(cro in 1:times_cross)
{ ######################
  #remove one from the population
  ########################
  L <- L_org
  X_l <- X_l_org
  X_u <- X_u_org
  X_lu <- rbind(X_l,X_u)
  lp  <- lp_org
  lu <- lu_org
  y_uun <- as.matrix(y_uun_org,ncol=1)
  labb <- which(y_uun==1)
  lagt <- which(y_uun==-1)
  one_l <- matrix(rep(1,times=n_l),ncol=1)
  one_u <- matrix(rep(1,times=n_u),ncol=1)
  one_lu <- matrix(rep(1,times=nn),ncol=1)
  I_l <- diag(n_l)
  I_u <- diag(n_u)
  I_beta<-diag(n_f)
  I_lu <- diag(nn)
  #initialize the objective variables

  


  KK_M11 <- lambda_2*I_beta + 2*C_2*t(X_u)%*%X_u + 2*t(X_lu)%*%L%*%X_lu + mu_1*t(X_l)%*%X_l
  KK_M12 <- 2*C_2*t(X_u)%*%one_u + 2*t(X_lu)%*%L%*%one_lu + mu_1*t(X_l)%*%one_l
  KK_M21 <- 2*C_2*t(one_u)%*%X_u + 2*t(one_lu)%*%L%*%X_lu + mu_1*t(one_l)%*%X_l 
  KK_M22 <- as.numeric(2*C_2*n_u+2*t(one_lu)%*%L%*%one_lu+mu_1*n_l)
  rbf <- rbfdot(sigma = 1)

  
  m<-ncol(X_lu)
  mm<-0.5*m*(m+1)+m
  X_lu_co<-fin_overall_convert(X_lu)
  X_lu_c<-re_format( X_lu_co,n_ccol=m)
  S_lu<-t(apply(matrix(1:(n_l+n_u),ncol=1),MARGIN=1,FUN=re_S,m=m,XX= X_lu_c,x= X_lu))
  X_test_co<-fin_overall_convert(TX2)
  X_test_c<-re_format( X_test_co,n_ccol=m)
  
  D_design<-fin_overall_convert(matrix(1,ncol=m,nrow=1))
  
  D_star_ini<-re_format( D_design,n_ccol=m)
  D_star_1<-design_final( D_star_ini)
  D_star<-design_format(D_star_1)
  D_star_0 <- matrix_0(D_star)
  
  G_re<-t(X_lu_c)%*%X_lu_c
  
  S_l<-as.matrix(S_lu[1:n_l,])
  S_u<-as.matrix(S_lu[(n_l+1):(n_l+n_u),])
  S_test<-t(apply(matrix(1:(nrow(TX2)),ncol=1),MARGIN=1,FUN=re_S,m=m,XX= X_test_c,x= TX2))
  
  
  p_alpha <- matrix(0.01,ncol=1,nrow=mm)
  p_beta_0<-0
  alpha_star <-matrix(data=p_alpha,ncol=1)
  beta_0<-p_beta_0
  


  h_ <- one_l -  (S_l %*% p_alpha + p_beta_0 * one_l)
  
  h_ <- matrix(data=h_,ncol=1)

  a_ <- one_u + (S_u %*% p_alpha + p_beta_0 * one_u)
  
  a_ <- matrix(data=a_,ncol=1)
  
  t_ <-  D_star%*%p_alpha
  
  f_ <- S_lu %*% p_alpha + p_beta_0 * one_lu
  
  v <- matrix(10,ncol=1,nrow=m*(m+1))
  
  u_h <- matrix(0,ncol=1,nrow=n_l)
  
  u_a <- matrix(0,ncol=1,nrow=n_u)
  
  q <- matrix(0,ncol=1,nrow=nn)
  
  
  
  

  #dual parameters
  u_h_0 <-rep(0.5,times=n_l)
  u_h <-matrix(data=u_h_0,ncol=1)


  RRecord<-matrix(nrow=0,ncol=nn+1+nrow(h_)+nrow(u_h)        )  #nn#
  sfk<-sample(labb,round(bv*(s_size-s_lp)))
  objj<-c()
  for(I_I in 1:iter_times)
  {  
    #####update alpha
     S_lu_star<- matrix_1(S_lu)
     S_l_star <-matrix_1(S_l)
     S_u_star <- matrix_1(S_u)
     I_m_0 <- matrix_0(diag(mm))
     G_re_0<- matrix_0(G_re)
     v_0<-vector_0(v)
     t_0<-vector_0(t_)
     M_111<-mu_3* t(D_star_0)%*% D_star_0+lambda_2*G_re_0+mu_1*t(S_l_star)%*%S_l_star+mu_2*t(S_u_star)%*%S_u_star+mu_4*t(S_lu_star)%*%S_lu_star
     M_222<-t(S_l_star)%*%u_h - t(S_u_star)%*%u_a -  t(D_star_0)%*%v_0 - t(S_lu_star)%*%q + mu_1*t(S_l_star)%*%(one_l-h_) - mu_2*t(S_u_star)%*%(one_u-a_) + mu_3* t(D_star_0)%*%t_0 + mu_4*t(S_lu_star)%*%f_
     M_solve<-solve(M_111)%*%M_222     
     alpha_star <- matrix( M_solve[1:mm],ncol=1)
     beta_0 <- as.numeric(M_solve[mm+1])
     
     
     #####update h
     update_lambda <- C_1/mu_1
     update_omega <- one_l+u_h/mu_1-(S_l%*% alpha_star +beta_0*one_l)
   
     h_ <- S_c(update_lambda,update_omega)
     
    
    #########update a
     a_<-1/(2*C_2+mu_2)*(u_a +  mu_2*(one_u+S_u%*% alpha_star+one_u*beta_0)   )
    
    
     
     ######### update t
     X_lasso<-diag(m*(m+1))
     Y_lasso<-D_star%*% alpha_star + v/mu_3
     i_for_t <- matrix(1:m,ncol=1)
     t_up<-apply(i_for_t,MARGIN=1,FUN=block_update,XX_bl=X_lasso,lambda=lambda_1,v=v,mu=mu_3,D_star=D_star,alpha_star=alpha_star)
     t_<-matrix(t_up,ncol=1,byrow=T)
     
     
      ######update of f
     f_<-solve(2*L + mu_4*diag(nn))%*%( q + mu_4*(S_lu%*%alpha_star+one_lu*beta_0) )
     
     #####dual
     
     u_h<-u_h+mu_1*( one_l -  (S_l %*%  alpha_star + beta_0 * one_l)-h_ )
     u_a<-u_a+mu_2*( one_u + (S_u %*%  alpha_star + beta_0 * one_u) -a_ )
     v <- v + mu_3*( D_star%*% alpha_star -t_ )
     q <- q + mu_4*( S_lu %*% alpha_star + beta_0 * one_lu - f_)
     
     
 
     
     
  

  }


  #######################
  #predict and error
  #######################
 # ssgss<-order(objj)[1]
  alpha_star<- cons_ver(vexc=alpha_star,fex=t_,desig=D_star)
  
  
  predict_test<-(sign(S_test %*% alpha_star + beta_0 ))
  error_test<-sum(abs(predict_test-TY2))/length(TY2)/2
  

  TP <- length(which(predict_test[TY2==1]==1))
  FP <- length(which(predict_test[TY2==-1]==1))
  TN <- length(which(predict_test[TY2==-1]==-1))
  FN <- length(which(predict_test[TY2==1]==-1))
  F1 <- 2*TP/(2*TP+FP+FN)    
  error_rate_overall<-c(error_rate_overall, F1) 
}
 

}



}}}}}


ecoli0.75<-error_rate_overall
print(ecoli0.75)
print(error_test)

