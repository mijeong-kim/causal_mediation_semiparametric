## jobs ###
library(devtools)
library(rootSolve)
library(mediation)
library(mvtnorm)
#install_url('https://cran.r-project.org/src/contrib/Archive/rmngb/rmngb_0.6-1.tar.gz')
library(rmngb)

### functions #############################
reg.m<-function(beta,x){ x%*%beta }
reg.m1<-function(beta,x){ t(x) }

################ Seff ################################################################
Seff<-function(theta,dta){ ## dta=list(x,y,v_cat)
  x<-as.data.frame(dta$x); y<-dta$y;
  x1<-x[,-1]; v_cat<-dta$v_cat; 

  p<-ncol(x1)
  for(i in 1:p){
    names(x1)[i]<-paste0("v",i)}
  
  nth<-length(theta);   nb<-nth-1
  beta<-theta[1:nb]; mu2<-theta[nth]
  e<-t(y-reg.m(beta,as.matrix(x))); n=length(e);
  e.dta<-data.frame(e=t(e),x1)
  cat<-which(v_cat==1); con<-which(v_cat==0)
  
  #form3<-as.formula(
  #  paste0("e^3~",paste0("s(v",con,")",collapse="+"),"+",paste0("v",cat,collapse="+"),collapse=""))
  #form4<-as.formula(
  #  paste0("e^4~",paste0("s(v",con,")",collapse="+"),"+",paste0("v",cat,collapse="+"),collapse=""))  
  #lm.e3<-gam(formula=form3,data=e.dta); 
  #lm.e4<-gam(formula=form4,data=e.dta); 
  
  form3<-as.formula(
    paste0("e^3~",paste0("v",con,collapse="^3+"),"^3+",
                  paste0("v",con,collapse="^2+"),"^2+",
                  paste0("v",1:p,collapse="+"),collapse="") 
    )
  
  form4<-as.formula(
    paste0("e^4~",paste0("v",con,collapse="^4+"),"^4+",
           paste0("v",con,collapse="^3+"),"^3+",
           paste0("v",con,collapse="^2+"),"^2+",
           paste0("v",1:p,collapse="+"),collapse="") 
  )
  
  lm.e3<-lm(formula=form3,data=e.dta); 
  lm.e4<-lm(formula=form4,data=e.dta); 
  mu3<-t(as.matrix(fitted(lm.e3))); mu4<-t(as.matrix(fitted(lm.e4)))
  t=e^2-mu2-mu3*e/mu2;  Et2=mu4-mu2^2-mu3^2/mu2;  mat1<-reg.m1(beta,x)
  
  Seff_beta<-NULL
  for(j in 1:nb){
    temp<-mat1[j,]*(e/mu2-mu3*t/mu2/Et2)
    Seff_beta<-rbind(Seff_beta,temp)}
  Seff_mu2<-t/Et2
  
  return(rbind(Seff_beta,Seff_mu2))
}


Seff_scaled<-function(theta,dta){
  rowSums(Seff(theta,dta))}

estcov<-function(theta,dta){
  x<-dta$x; y<-dta$y
  n=length(y); p=length(theta)
  f=Seff(theta,dta)
  B=f%*%t(f)/n
  A=matrix(0,nrow=p,ncol=p)
  delta=1e-3
  for(i in 1:p){
    d=rep(0,p)
    d[i]=theta[i]*delta
    A[,i]=rowMeans(Seff(theta+d,dta)-Seff(theta-d,dta))/d[i]/2
  }
  M=solve(A)%*%B%*%t(solve(A))
  return(M)
}

##########################
Seff1<-function(theta,dta){
  x<-as.data.frame(dta$x); y<-dta$y;
  nth<-length(theta);
  nb<-nth-1
  beta<-theta[1:nb]; mu2<-theta[nth]
  e<-t(y-reg.m(beta,as.matrix(x)))
  n=length(e)
  h=1.06*sd(e)*n^(-1/5);
  fhat=rep(0,n)
  fdhat=rep(0,n)
  for(i in 1:n){
    fhat[i]=sum(dnorm((e[i]-e)/h))/(n*h);
    fdhat[i]=sum(-(e[i]-e)*dnorm((e[i]-e)/h))/(n*h^3);
  }
  mu3<-sum(e^3)/n;
  mu4=sum(e^4)/n;
  t=e^2-mu2-mu3*e/mu2;
  Et2=mu4-mu2^2-mu3^2/mu2;
  f1f=fdhat/fhat;
  
  mat1<-(reg.m1(beta,x)-rowMeans(reg.m1(beta,x)))
  mat2<-rowMeans(reg.m1(beta,x))
  
  Seff_beta<-NULL
  for(j in 1:nb){
    temp<--f1f*mat1[j,]+(e/mu2-mu3*t/mu2/Et2)*mat2[j]
    Seff_beta<-rbind(Seff_beta,temp)
  }
  Seff_mu2<-t/Et2
  rbind(Seff_beta,Seff_mu2)
}

Seff_scaled1<-function(theta,dta){
  rowSums(Seff1(theta,dta))}

estcov1<-function(theta,dta){
  x<-dta$x; y<-dta$y
  n=length(y)
  p=length(theta)
  f=Seff1(theta,dta)
  B=f%*%t(f)/n
  A=matrix(0,nrow=p,ncol=p)
  delta=1e-3
  for(i in 1:p){
    d=rep(0,p)
    d[i]=theta[i]*delta
    A[,i]=rowMeans(Seff1(theta+d,dta)-Seff1(theta-d,dta))/d[i]/2
  }
  M=solve(A)%*%B%*%t(solve(A))
  return(M)
}

##########################
##########################

causal.eff<-function(J=100,K=1000,  med.para, out.para, med.data, out.data){
  ACME1<-rep(NA,J); ACME0<-rep(NA,J);   ACME<-rep(NA,J);
  ADE1<-rep(NA,J); ADE0<-rep(NA,J);   ADE<-rep(NA,J);
  TE<-rep(NA,J);

  mu.med<-med.para$th;  Sigma.med<-med.para$cov
  y.med<-med.data$y; n=length(y.med)
  x.med.mat<-as.matrix(data.frame(ones=rep(1,n),med.data[,-1]));   
  th.med.J<-rmvnorm(J,mu.med,Sigma.med)
  e.med<-y.med-x.med.mat%*%as.vector(mu.med)

  mu.out<-out.para$th;  Sigma.out<-out.para$cov
  y.out<-out.data$y; 
  x.out.mat<-as.matrix(data.frame(ones=rep(1,n),out.data[,-1]));   
  th.out.J<-rmvnorm(J,mu.out,Sigma.out)
  e.out<-y.out-x.out.mat%*%as.vector(mu.out)
  
  x.med.t0.dta<-x.med.mat;  x.med.t0.dta[,"treat"]<-0; x.med.t0<-as.matrix(x.med.t0.dta)
  x.med.t1.dta<-x.med.mat;  x.med.t1.dta[,"treat"]<-1; x.med.t1<-as.matrix(x.med.t1.dta)
  x.out.t1.M1.dta<-x.out.mat; x.out.t1.M1.dta[,"treat"]<-1;   x.out.t1.M0.dta<- x.out.t1.M1.dta
  x.out.t0.M1.dta<-x.out.mat; x.out.t0.M1.dta[,"treat"]<-0;   x.out.t0.M0.dta<- x.out.t0.M1.dta
  
  M.t0<-matrix(NA,nrow=n,ncol=K);   M.t1<-matrix(NA,nrow=n,ncol=K)
  Y.t1.M1<-matrix(NA,nrow=n,ncol=K);  Y.t1.M0<-matrix(NA,nrow=n,ncol=K)
  Y.t0.M1<-matrix(NA,nrow=n,ncol=K);  Y.t0.M0<-matrix(NA,nrow=n,ncol=K)
  
  for(j in 1:J){
    for(k in 1:K){
      M.t0[,k]<-x.med.t0%*%as.vector(th.med.J[j,])+rDist(e.med)(1)
      M.t1[,k]<-x.med.t1%*%as.vector(th.med.J[j,])+rDist(e.med)(1)
      ##
      x.out.t1.M1.dta[,"M"]<-M.t1[,k];  
      x.out.t1.M1<-as.matrix(x.out.t1.M1.dta)
      Y.t1.M1[,k]<-x.out.t1.M1%*%as.vector(th.out.J[j,])+rDist(e.out)(1)
      ##
      x.out.t1.M0.dta[,"M"]<-M.t0[,k];  
      x.out.t1.M0<-as.matrix(x.out.t1.M0.dta)
      Y.t1.M0[,k]<-x.out.t1.M0%*%as.vector(th.out.J[j,])+rDist(e.out)(1)        
      ##
      x.out.t0.M1.dta[,"M"]<-M.t1[,k];  
      x.out.t0.M1<-as.matrix(x.out.t0.M1.dta)
      Y.t0.M1[,k]<-x.out.t0.M1%*%as.vector(th.out.J[j,])+rDist(e.out)(1)
      ##
      x.out.t0.M0.dta[,"M"]<-M.t0[,k];  
      x.out.t0.M0<-as.matrix(x.out.t0.M0.dta)
      Y.t0.M0[,k]<-x.out.t0.M0%*%as.vector(th.out.J[j,])+rDist(e.out)(1)  
    }
    ACME1[j]<-sum(Y.t1.M1-Y.t1.M0)/(n*K); ACME0[j]<-sum(Y.t0.M1-Y.t0.M0)/(n*K) 
    ACME[j]<-(ACME1[j]+ACME0[j])/2;
    ADE1[j]<-sum(Y.t1.M1-Y.t0.M1)/(n*K);  ADE0[j]<-sum(Y.t1.M0-Y.t0.M0)/(n*K)
    ADE[j]<-(ADE1[j]+ADE0[j])/2;
  }
  TE<-ACME+ADE;
  delta<-quantile(ACME, probs=c(0.025,0.5,0.975))
  zeta<-quantile(ADE, probs=c(0.025,0.5,0.975))
  te<-quantile(TE, probs=c(0.025,0.5,0.975))
  m<-c(mean(ACME),mean(ADE),mean(te))
  list(ACME=delta,ADE=zeta,Total=te,Estimate=m)
}

##########################
##########################
#### plot for ACME, ADE, Total Effect ####
plot_effects<-function(causal.eff,main_text=""){
  effect <- c("Total\n Effect","ADE","ACME")
  estimate = causal.eff$Estimate[3:1] # 각 효과의 추정값
  lower = c(causal.eff$Total[1], causal.eff$ADE[1], causal.eff$ACME[1])    # 신뢰구간 하한
  upper = c(causal.eff$Total[3], causal.eff$ADE[3], causal.eff$ACME[3])
  
  # 빈 플롯 생성
  #plot(x = NULL, y = NULL, xlim = c(min(lower), max(upper)), ylim = c(0.5, length(effect)+0.5),
  plot(x = NULL, y = NULL, xlim=c(-0.19,0.09), ylim = c(0.5, length(effect)+0.5),
       xlab = "", ylab = "", axes = FALSE, main=main_text)
  
  
  # y축 라벨 추가
  axis(2, at = 1:length(effect), labels = effect, las = 1)
  
  # x축 라벨 추가
  axis(1)
  
  # 참조선을 그리기 위해 x=0에 수직선 추가
  abline(v = 0, lty = 2 )
  
  # 각 효과의 추정값을 점으로 표시
  points(estimate, 1:length(effect), pch = 19, cex=0.8)
  
  # 신뢰구간을 선으로 표시
  for (i in 1:length(effect)) {
    segments(x0 = lower[i], x1 = upper[i], y0 = i, y1 = i , lwd=1.5)
  }
  
  # 플롯 테두리 추가 (선택 사항)
  box()
}
