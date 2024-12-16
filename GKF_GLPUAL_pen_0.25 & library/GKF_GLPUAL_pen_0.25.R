#####################
# environment
#####################
#packages
library(mvtnorm)
library(expm)
#times of iteration for model training
iter_times=1100
#tool functions
source('F:/ucl/pr/result/script/EX_fun_PU.R')

###################################
#Load original dataset
###################################
data_1<-read.csv("F:/ucl/pr/dataset/112/ngood/pendigits.csv",header = F)
data_class<-data_1[,17]
data_original<-rbind(data_1[data_class==1,][1:200,],data_1[data_class==8,][1:200,],data_1[data_class==4,][1:400,])

###################################
#Random Seeds
###################################
ran_start_1<-12235
ran_start_2<-12207
set.seed(ran_start_1)


###############################################
#The setting of hyper-parameters
###############################################
#hyper-parameters in the objective function
lambda_1 = 0.1
lambda_2 = 0.17
C_1 = 1
C_2 = 0.01
mu_1 = 0.01
#Augmented Lagrangian coefficients
mu_2<-0.01
mu_3<-0.01
mu_4<-0.01
mu_3_star<-0.001

##########################################
#Add irrelevant attributes to dataset.
##########################################
mg_meaningless<-rnorm(nrow(data_original), mean =0, sd = 1)
for(mgm in 1:40)
  {
  mg_m1<-rnorm(nrow(data_original), mean =0, sd = 1)
  mg_meaningless<-cbind(mg_meaningless,mg_m1)
  }
data_original<-cbind(data_original,mg_meaningless)
PP_cro<-c() 
F1_overall <- c()


set.seed(ran_start_2)
###################################
#Construction of PU data
###################################
ranran<-sample(1:nrow(data_original))
data1<-data_original[ranran,]
data1<- data1[sample(1:nrow(data1)),]
data11<-data1[,c(1:16,18:33)]
data12<-data1[,17]

#Define the positive class and negative class
data12[data12==1] <- 1
data12[data12==8] <- 1
data12[data12==4] <- -1
data12<-as.numeric(data12)
fre_label=0.75
p_trte=0.7
re_sig=0.001
p.pi<-length(data12[data12==1])/length(data12)
  
#Remote labels to construct PU data
data13 <- data12
ta<-sample(1:length(data13[data13==1] ))[1:(fre_label*length(data13[data13==1]))]
data12[data13==1][ta]<--1
data11<-data.matrix(data11, rownames.force = NA)
data11<-data11
data112<-data11
data_non<-data11[data12==-1,]
data_fra<-data11[data12==1,]
data_fra<-data_fra


#Split between training and test set.
Trainingset_1<-data13[1:(p_trte*length(data13))]
Training_Y<-matrix(Trainingset_1,ncol=1)
n_f<-ncol(data11)
Trainingset_2<-data11[1:(p_trte*length(data13)),]
Training_X<-matrix(Trainingset_2,ncol=n_f)
Testset<-c(data13[(p_trte*length(data13)+1):length(data13)])
Test_Y<-matrix(Testset,ncol=1)
Testset4<-rbind(data11[(p_trte*length(data13)+1):length(data13),])
Test_X<-matrix(Testset4,ncol=n_f)

#sample size, positive and negative
s_size <- length(which(Training_Y==1))
s_size_n <-length(which(Training_Y==-1))
#size of labelled positive training set
s_lp <-length(which(data12[1:(0.7*length(data13))]==1))

#labelled instances in positives and unlabelled instances
lp <- which(data12[1:(0.7*length(data13))]==1)
lu <-intersect( which((Training_Y==1)), which(data12[1:(0.7*length(data13))]==-1))

#Positive matrix
X_P <- Training_X[Training_Y==1,]
#Negative matrix
X_N <- Training_X[Training_Y=-1,]
#labelled and unlabelled matrices 
X_l <-  Training_X[data12[1:(0.7*length(data13))]==1,]
X_u <- Training_X[data12[1:(0.7*length(data13))]==-1,]
#matrices for Y
Y_uun1 <- matrix(1,ncol=1,nrow=s_size-s_lp)
Y_uun2 <- matrix(-1,ncol=1,nrow=s_size_n)
y_uun<- rbind(Y_uun1,Y_uun2)


#size of positive and negative
n_l <- nrow(X_l)
n_u <- nrow(X_u)
#overall sample size
nn <- n_l+n_u
X_lu<-rbind(X_l,X_u)
#constant vectors and matrices
one_l <- matrix(rep(1,times=n_l),ncol=1)
one_u <- matrix(rep(1,times=n_u),ncol=1)
one_lu <- matrix(rep(1,times=nn),ncol=1)
I_l <- diag(n_l)
I_u <- diag(n_u)
I_beta<-diag(n_f)
I_lu <- diag(nn)


###########################
#the similar matrix L
###########################
#find neighbors
theta_K<-1
W<-matrix(0,nrow=nn,ncol=nn)
Phi_lu<-matrix(0,nrow=nn,ncol=nn)
#calculate the similar matrix L
for(j in 1:nn)
{ 
    for(i in 1:nn)
    {
       W[i,j]<-exp(- norm( matrix(X_lu[i,]-X_lu[j,],ncol = 1) ,type="2")^2 /2/theta_K)
    }
}
diag(W)<-rep(0,times=nn)
D<-diag(rowSums(W))
L <- D-W

###############################
#Labels matrices in the model
###############################
y_l <- matrix(rep(1,times=n_l),ncol=1)
y_u <- matrix(rep(-1,times=n_u),ncol=1)
y_lu <- rbind(y_l,y_u)
Y_l <- matrix(0,ncol=n_l,nrow=n_l)
Y_u <- matrix(0,ncol=n_u,nrow=n_u)
Y_lu <- matrix(0,ncol=nn,nrow=nn)
diag(Y_l) <- y_l
diag(Y_u) <- y_u
diag(Y_lu) <- y_lu
  
##########################################################
#Change the form of the model matching to the linear case
##########################################################
m<-ncol(X_lu)
mm<-0.5*m*(m+1)+m
X_lu_co<-fin_overall_convert(X_lu)
X_lu_c<-re_format( X_lu_co,n_ccol=m)
S_lu<-t(apply(matrix(1:(n_l+n_u),ncol=1),MARGIN=1,FUN=re_S,m=m,XX= X_lu_c,x= X_lu))
X_test_co<-fin_overall_convert(Test_X)
X_test_c<-re_format( X_test_co,n_ccol=m)
D_design<-fin_overall_convert(matrix(1,ncol=m,nrow=1))
D_star_ini<-re_format( D_design,n_ccol=m)
D_star_1<-design_final( D_star_ini)
D_star<-design_format(D_star_1)
D_star_0 <- matrix_0(D_star)
G_re<-t(X_lu_c)%*%X_lu_c
S_l<-as.matrix(S_lu[1:n_l,])
S_u<-as.matrix(S_lu[(n_l+1):(n_l+n_u),])
S_test<-t(apply(matrix(1:(nrow(Test_X)),ncol=1),MARGIN=1,FUN=re_S,m=m,XX= X_test_c,x= Test_X))

###################################################   
#initialize the model parameters and dual variables
###################################################
#model parameters
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

#dual variables
u_h <- matrix(0,ncol=1,nrow=n_l)
u_a <- matrix(0,ncol=1,nrow=n_u)
q <- matrix(0,ncol=1,nrow=nn)
u_h_0 <-rep(0.5,times=n_l)
u_h <-matrix(data=u_h_0,ncol=1)
##################################
#Model Training
##################################
for(I_I in 1:iter_times)
    { 
    ##############################
    #Update alpha
    ###############################
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
   
    ########################################
    #Update of h
    ########################################
    update_lambda <- C_1/mu_1
    update_omega <- one_l+u_h/mu_1-(S_l%*% alpha_star +beta_0*one_l)
    h_ <- S_c(update_lambda,update_omega)
   
    #########################################    
    #Update of a
    ########################################
    a_<-1/(2*C_2+mu_2)*(u_a +  mu_2*(one_u+S_u%*% alpha_star+one_u*beta_0)   )
  
  
    ########################################
    # Update of t
    ########################################
    X_lasso<-diag(m*(m+1))
    Y_lasso<-D_star%*% alpha_star + v/mu_3
    i_for_t <- matrix(1:m,ncol=1)
    t_up<-apply(i_for_t,MARGIN=1,FUN=block_update,XX_bl=X_lasso,lambda=lambda_1,v=v,mu=mu_3,D_star=D_star,alpha_star=alpha_star)
    t_<-matrix(t_up,ncol=1,byrow=T)
   
    ########################################
    #Update of f
    ########################################
    f_<-solve(2*L + mu_4*diag(nn))%*%( q + mu_4*(S_lu%*%alpha_star+one_lu*beta_0) )
    
    ########################################
    #Update of dual variables
    ########################################
    u_h<-u_h+mu_1*( one_l -  (S_l %*%  alpha_star + beta_0 * one_l)-h_ )
    u_a<-u_a+mu_2*( one_u + (S_u %*%  alpha_star + beta_0 * one_u) -a_ )
    v <- v + mu_3*( D_star%*% alpha_star -t_ )
    q <- q + mu_4*( S_lu %*% alpha_star + beta_0 * one_lu - f_)
    }
        
#############################################################  
#Change the form of model parameters back to the initial one
#############################################################
alpha_star<- cons_ver(vexc=alpha_star,fex=t_,desig=D_star)

######################################################
#Predict the labels for the test set
####################################################
predict_test<-sign(S_test %*% alpha_star + beta_0 )

####################################################
#Calculate F1 Score adn error rate on the test set
####################################################
TP <- length(which(predict_test[Test_Y==1]==1))
FP <- length(which(predict_test[Test_Y==-1]==1))
TN <- length(which(predict_test[Test_Y==-1]==-1))
FN <- length(which(predict_test[Test_Y==1]==-1))
F1 <- 2*TP/(2*TP+FP+FN)    
F1_overall<-c(F1_overall, F1)
error_test<-sum(abs(predict_test-Test_Y))/length(Test_Y)/2

###################################
#Output F1 Score on the test set
###################################
print("F1 Score")
print(F1_overall)
print("error rate")
print(error_test)



