if(iter_times%%10!=0)
{h_00<-h_
a_00<-a_
t_00<-t_
f_00<-f_}

eps<-0.001
eps2<-0.25

norm_1<-norm(t(cbind(S_l,one_l))%*%u_h,type="2")
norm_2<-norm(t(cbind(S_u,one_u))%*%u_a,type="2")
norm_3<-norm(t(D_star)%*%v,type="2")
norm_4<-norm(t(cbind(S_lu,one_lu))%*%q,type="2")

thre_stop_1<-eps2*sqrt(mm+1)+eps*norm_1
thre_stop_2<-eps2*sqrt(mm+1)+eps*norm_2
thre_stop_3<-eps2*sqrt(mm+1)+eps*norm_3
thre_stop_4<-eps2*sqrt(mm+1)+eps*norm_4

vari_1<-norm(mu_1*t(cbind(S_l,one_l))%*%(h_-h_00))
vari_2<-norm(mu_2*t(cbind(S_u,one_u))%*%(a_-a_00))
vari_3<-norm(mu_3*t(D_star)%*%(t_-t_00),type="2")
vari_4<-norm(mu_4*t(cbind(S_lu,one_lu))%*%(f_-f_00))

if((vari_1<thre_stop_1)&(vari_2<thre_stop_2)&(vari_3<thre_stop_3)&(vari_4<thre_stop_4)){print("converged")}