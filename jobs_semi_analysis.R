## jobs ###
library(devtools)
library(rootSolve)
library(mediation)
library(mvtnorm)
library(gam)
#install_url('https://cran.r-project.org/src/contrib/Archive/rmngb/rmngb_0.6-1.tar.gz')
library(rmngb)


################################################################################
################################################################################
################################################################################
### Calculation
### DATA : jobs
data(jobs)
med.lm <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
out.lm <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs) ### No interaction
# Estimation via quasi-Bayesian approximation
contcont <- mediate(med.lm, out.lm, sims=50, treat="treat", mediator="job_seek")
summary(contcont)
plot(contcont)

# Estimation via nonparametric bootstrap
#contcont.boot <- mediate(med.lm, out.lm, boot=TRUE, sims=50, treat="treat", mediator="job_seek")
#summary(contcont.boot)
#plot(contcont.boot)

res.med<-residuals(med.lm); shapiro.test(res.med)
res.out<-residuals(out.lm); shapiro.test(res.out)

ACME.lm<-med.lm$coefficients["treat"]*out.lm$coefficients["job_seek"]
ADE.lm<-out.lm$coefficients["treat"]
TE.lm<-ACME.lm+ADE.lm
c(ACME.lm, ADE.lm, TE.lm)

##################################
####### equation 1################
eq1.lm <- lm(depress2 ~ treat + econ_hard + sex + age, data=jobs)

################################################################################
################################################################################
####### semiparametric #####
### Mediation model ###
y.med<-jobs$job_seek
n=length(y.med)
x.med.dta<-data.frame(ones=rep(1,n),treat=jobs$treat, 
                      econ_hard=jobs$econ_hard, 
                      sex=jobs$sex, 
                      age=jobs$age)
x.med<-as.matrix(x.med.dta)
v.cat<-c(1,0,1,0)
dta.med<-list(x=x.med,y=y.med,v_cat=v.cat)

#### mediation, E(e|X,t) ###
#alpha=0.93; 0.5264107 0.4447724 
#alpha=0.91; 0.5458160 0.5666928 
#alpha=0.89; 0.5752367 0.5776999
#alpha=0.83; 0.5267161 0.5620824 ; p-value 유의하지 않음
#alpha=0.81; 0.5264105 0.5339081 p-value 유의하지 않음
#alpha=0.74; 1.096664 1.063715; 모두 유의함
#alpha=0.67; 0.5660368 0.5906906; 모두 유의함
#alpha=0.60; 0.5548939 0.5139363; 전반적으로 유의함


e.med<-y.med-x.med%*%coef(med.lm)
th0.med<-c(med.lm$coefficients,var(e.med)*0.91)
result.med<-multiroot(Seff_scaled,start=th0.med, dta=dta.med)
th.semi.med<-result.med$root  
nth<-length(th.semi.med); nb=nth-1;
res.semi.med<-y.med-x.med%*%th.semi.med[1:nb]
hist(res.semi.med)

M.med=estcov(th.semi.med,dta.med)/n
format(M.med[1:5,1:5],scientific=FALSE)
se.med<-sqrt(diag(M.med))

z.ratio<-th.semi.med/se.med; 
pz<-(1-pnorm(abs(z.ratio)))*2; pz1<-format(pz,scientific=FALSE)
list(th=th.semi.med, se=se.med, z.ratio=z.ratio, p.value=pz1)

### Comparison of variance of residuals
c(var(res.semi.med),  th.semi.med[6] )



### outcome model (no interaction)
y.out<-jobs$depress2
x.out.dta<-data.frame(ones=rep(1,n),treat=jobs$treat, job_seek=jobs$job_seek, 
                      econ_hard=jobs$econ_hard, sex=jobs$sex, age=jobs$age)
x.out<-as.matrix(x.out.dta)
dta.out<-list(x=x.out,y=y.out, v_cat=c(1,0,0,1,0))
e.out<-y.out-x.out%*%coef(out.lm)

##
#alpha=0.98; 0.4743304 0.4503614; p-value 모두 0
#alpha=0.94; 0.7907020 0.6931233;p-value 모두 0
#alpha=0.93; 0.4427864 0.4700978;p-value 모두 0
#alpha=0.90; 0.5977470 0.6294083 ; p-value 모두 0
#alpha=0.88; 0.6317042 0.6306601 ; p-value 모두 0
#alpha=0.85; 0.4048038 0.4381204
#alpha=0.83;  0.7039732 0.6419791  ; p-value 모두 0
#alpha=0.82; 0.3801126 0.4546611  ; p-value 모두 0
#alpha=0.80; 0.3822170 0.3115497; p-value 모두 0
#alpha=0.77; 0.3761801 0.412104 good!!!!




th0.out<-c(out.lm$coefficients,var(e.out)*0.77)
result.out<-multiroot(Seff_scaled,start=th0.out, dta=dta.out, maxiter = 200)
th.semi.out<-result.out$root  
nth<-length(th.semi.out); nb=nth-1;
res.semi.out<-y.out-x.out%*%th.semi.out[1:nb]
hist(res.semi.out)

M.out=estcov(th.semi.out,dta.out)/n
format(M.out[1:6,1:6],scientific=FALSE)
se.out<-sqrt(diag(M.out)); 
c(var(res.semi.out),  th.semi.out[7] )

z.ratio<-th.semi.out/se.out
pz<-(1-pnorm(abs(z.ratio)))*2;  pz1<-format(pz,scientific=FALSE);
list(th=th.semi.out, se=se.out, z.ratio=z.ratio, p.value=pz1)

format(th.semi.out,scientific=FALSE)





############################################################################################
############################################################################################
####### semiparametric #### e,X :indep ########################
### Mediation model ###
dta1.med<-list(x=x.med,y=y.med)
e.med<-y.med-x.med%*%coef(med.lm)
th0.med<-c(med.lm$coefficients,var(e.med))
result1.med<-multiroot(Seff_scaled1,start=th0.med, dta=dta1.med)
th1.semi.med<-result1.med$root  
nth<-length(th1.semi.med); nb=nth-1;
res1.semi.med<-y.med-x.med%*%th1.semi.med[1:nb]
hist(res.semi.med)

### Numerical method for calculation of SD
M1.med=estcov1(th1.semi.med,dta1.med)/n
format(M1.med[1:5,1:5],scientific=FALSE)
se1.med<-sqrt(diag(M1.med))

z.ratio<-th1.semi.med/se1.med
pz<-(1-pnorm(abs(z.ratio)))*2;  pz1<-format(pz,scientific=FALSE);
list(th=th1.semi.med, se=se1.med, z.ratio=z.ratio, p.value=pz1)


### Comparison of variance of residuals
c(var(res1.semi.med),  th1.semi.med[6] )


### outcome model (no interaction)
dta1.out<-list(x=x.out,y=y.out)
e.out<-y.out-x.out%*%coef(out.lm)
th0.out<-c(out.lm$coefficients,var(e.out))
result1.out<-multiroot(Seff_scaled1,start=th0.out,dta=dta1.out)
th1.semi.out<-result1.out$root  
nth<-length(th1.semi.out); nb=nth-1;
res1.semi.out<-y.out-x.out%*%th1.semi.out[1:nb]
hist(res1.semi.out)

M1.out=estcov1(th1.semi.out,dta1.out)/n
format(M1.out[1:6,1:6],scientific=FALSE)
se1.out<-sqrt(diag(M1.out))

z.ratio<-th1.semi.out/se1.out
pz<-(1-pnorm(abs(z.ratio)))*2;  pz1<-format(pz,scientific=FALSE);
list(th=th1.semi.out, se=se1.out, z.ratio=z.ratio, p.value=pz1)



### Causal effect ###
##causal.eff(J=100,K=1000,  med.para, out.para, med.data, out.data)
#  mu.med<-med.para$th;  Sigma.med<-med.para$cov
med.para<-list(th=th.semi.med[1:5],cov=M.med[1:5,1:5])
out.para<-list(th=th.semi.out[1:6],cov=M.out[1:6,1:6])
med.dta<-data.frame(y=y.med,x.med.dta[,-1])
out.dta<-data.frame(y=y.out,x.out.dta[,-1])
names(out.dta)[3]<-"M"


### Causal effect ###
##causal.eff(J=100,K=1000,  med.para, out.para, med.data, out.data)
#  mu.med<-med.para$th;  Sigma.med<-med.para$cov
med1.para<-list(th=th1.semi.med[1:5],cov=M1.med[1:5,1:5])
out1.para<-list(th=th1.semi.out[1:6],cov=M1.out[1:6,1:6])
#med.dta<-data.frame(y=y.med,x.med.dta[,-1])
#out.dta<-data.frame(y=y.out,x.out.dta[,-1])
#names(out.dta)[3]<-"M"



#set.seed(62)
#set.seed(25)
set.seed(95)
contcont <- mediate(med.lm, out.lm, sims=100, treat="treat", mediator="job_seek")
semi.c_eff<-causal.eff(J=100,K=1000,  med.para, out.para, med.dta, out.dta)
semi1.c_eff<-causal.eff(J=100,K=1000,  med1.para, out1.para, med.dta, out.dta)

summary(contcont)

######  plot 그리기 ####
setwd("C:\\Users\\user\\Dropbox\\2_research\\0_causal_semipara\\Latex\\figure")
pdf("causal_effects_output.pdf", width=8.27, height=11.69)
par(mfrow=c(3,1))
plot(contcont,xlim=c(-0.19,0.09),main="(a) OLS")
plot_effects(semi.c_eff,
             expression(bold("(b) Semiparametric: ") * "no assumption for " * epsilon))
plot_effects(semi1.c_eff,
             expression(bold("(c) Semiparametric: ") * epsilon
                             * " independent with (" * italic(X) * "," * italic(T) * ")"))

dev.off()
