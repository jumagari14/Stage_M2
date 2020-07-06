##### fit de ks et kd #######################################
# on dispose de rna(t) et de n mesures (ti,pi)
# dp/dt=ks*rna-(kd+mu)*p
# p(0)=p0
# ###########################################################
setwd("../data_kiwi/")
library(deSolve)
library(rlang)# utilisation de is empty
library(pracma)
library(minpack.lm)

############## avec mu ##################################
# donn?es poids
poids <- read.csv("poids_kiwi.csv",header = T,col.names = c("t","w"))
attach(poids)
plot(t,w,main="poids tomate")
# fit du poids par double sigmoide
sigmoide=function(t,par){
  s=par[4]+par[1]/(1+exp(-par[2]*(t-par[3]))) +par[5]/(1+exp(-par[6]*(t-par[7])))
  return(s)}
#

initialisation=c(a=48.4,b=0.27,c=19.62,d=0.1,e=80.27,f=0.11,g=47.82)
curve(sigmoide(x,initialisation),0,55,ylim=c(0,150))
points(t,w,col="red")
model_sig=nls(w~d+a/(1+exp(-b*(t-c))) +e/(1+exp(-f*(t-g))),start = initialisation,algorithm = "port")
# fit du poids par double sigmoide
mu=function(t,par){
  s=par[4]+par[1]/(1+exp(-par[2]*(t-par[3]))) +par[5]/(1+exp(-par[6]*(t-par[7])))
  ds=par[1]*par[2]*exp(-par[2]*(t-par[3]))/(1+exp(-par[2]*(t-par[3])))^2
  par[5]*par[6]*exp(-par[6]*(t-par[7]))/(1+exp(-par[6]*(t-par[7])))^2
  mu=ds
  return(mu)}
#
valp=coef(model_sig)
curve(mu(x,valp),0,50,ylim=c(0,6))
##############################################################
# choix rna courbe
rna<-function(t) {
  #out=-10/55^2*(t-55/2)^2+2.5
  out=2*exp(-(t-55/2)^2/500)
  return(out)}

rna<-function(t) {
  #out=-10/55^2*(t-55/2)^2+2.5
  out=2*exp(-t/5)-1/70*t+11/14
  return(out)}
rna<-function(t) {
  #out=-10/55^2*(t-55/2)^2+2.5
  #out=0.003124*t^2-0.217440*t+4.009098
  out=0.003559*t^2-0.233527*t+4.078663
  return(out)}
dt=7
t1=seq(0,55,dt)
ti=sort(c(t1,t1,t1))
test_rna<-rna(ti)
test_rna<-test_rna/mean(test_rna)
coef=polyfit(ti,log(test_rna),3)
# g?n?ration des donn?es par r?solution ode


plot(ti,test_rna,"l",col="green")
lines(ti,polyval(coef,t),"l",col="red")
# cr?ation donn?es avec prise en compte de mu
##############################################################
fmu=function(t,p,par){
  dp=par["ks"]*exp(polyval(coef,t))-(par["kd"]+mu(t,valp))*p
  return(list(dp))}
p0=3
ks_vrai=0.5
kd_vrai=0.2
parInit<-list("start_prot"=p0,"ks"=ks_vrai,"kd"=kd_vrai)
parInit2<-list("start_prot"=p0,"ks"=ks_vrai,"kd"=kd_vrai*0.1)
parInit3<-list("start_prot"=p0,"ks"=ks_vrai*5,"kd"=kd_vrai*5)
parInit4<-list("start_prot"=p0,"ks"=ks_vrai*10,"kd"=kd_vrai*10)
listPar<-list("par1"=parInit,"par2"=parInit2,"par3"=parInit3,"par4"=parInit4)
# r?solution EDO pour ti=0,7,14,..
# ajout de bruit de 20%
bruit=1

Noisify <- function(data,bruit) {
  #  noise <- runif(length(data), -bruit, bruit)
  noise<-rnorm(length(data),0,bruit)
  noisified <- data+noise
  noisified[noisified<0]=abs(data[noisified<0])
  # si la valeur bruit?e est n?gative, on diminue le bruit 
  #noisified[noisified<0]=data[noisified<0]+runif(length(noisified[noisified<0]), -data[noisified<0], data[noisified<0])
  #noisified[noisified<0]=data[noisified<0]+rnorm(length(noisified[noisified<0]), 0,data[noisified<0] )
  return(noisified)
}
pdf("Fake data test.pdf")
res_list<-lapply(listPar, function(el){
  tryCatch({
  dt=7
  t1=seq(0,55,dt)
  ti=sort(c(t1,t1,t1))
  sol=ode(y=p0,times=ti,func=fmu,parms=c(ks=el[["ks"]],kd=el[["kd"]]))
  
  solb<-matrix(nrow=length(sol[,2]),ncol=2)
  solb[,1]=sol[,1]
  bruit=min( sqrt(mean(sol[,2])/6),0.5)
  val_prot<-Noisify(sol[,2],bruit)
  solb[,2]=val_prot
  el[["Noisy sol"]]=solb
  el[["Original sol"]]=solb
  # resolution pour ti=0,7,14,...
  #sol=ode(y=p0,times=t,func=f,parms=c(ks=0.31,kd=0.17))
  plot(solb,ylim=c(0,max(sol[,2])),col="green")
  points(sol,ylim=c(0,max(sol[,2])),col="blue")
  curve(rna(x),0,55,col="blue",add=TRUE)
  legend("topleft",c("pts bruit?s","points exacts","rna"),lty=c(NA,NA,1),pch=c(1,1,NA),col =c("green","blue","blue"),cex=0.5)
  solb=data.frame(solb)
  names(solb) <- c("t", "p")
  attach(solb)
  
  
  ##################################################################"
  # nls pour retrouver ks et kd
  ##################################################################
  
  #nls.control(maxiter = 1000, tol = 1e-03, minFactor = 1/1024,
  #printEval = FALSE, warnOnly = T)
  #ksmin=ks_vrai*0.25#
  ksmin=0
  ksmax=el$ks*40
  #kdmin=kd_vrai*0.25
  kdmin=0
  kdmax=el$kd*40
  lb=c(ksmin,kdmin,0);ub=c(ksmax,kdmax,2*p0)
  # starting value de kd: max(...) to avoid negative value
  kdinit=max(log(solb[1,2]/solb[4,2])/7,0.001)
  # starting value ks
  ksinit=max(solb[4,2]/solb[4,1]*kdinit,0.001)
  initialisation=list(ks=ksinit,kd=kdinit,p00=p0)
  fit<-nls(solb$p~ode(y=p00,times=sol[,1],func=fmu,parms=c(ks=ks,kd=kd))[,2],data=list(solb),start=initialisation,algorithm="port",control=list(minFactor=1e-5,maxiter=1000),lower=lb,upper=ub)
  # fit2<-nlsLM(solb$p~ode(y=p00,times=sol[,1],func=fmu,parms=c(ks=ks,kd=kd)),data=list(solb),start=initialisation,control=list(minFactor=1e-5,maxiter=1000),lower=lb,upper=ub)
  summary(fit)
  # summary(fit2)
  #fit<-nls(solb$p~ode(y=coef(fit)[3],times=sol[,1],func=f,parms=c(ks=ks,kd=kd))[,2],data=list(solb),start=coef(fit),algorithm="port",control=list(minFactor=1e-5,maxiter=1000),lower=lb,upper=ub)
  el[["coefs"]]=coef(fit)
  theo<-ode(y=coef(fit)[3],times=sol[,1],func=fmu,parms=coef(fit))
  hauteur=max(max(sol[,2]),max(solb[,2]))
  plot(theo[,1],theo[,2],type="l",main="Fit Test ", col="green",ylim=c(0,14),xlab="t (days)", ylab="Proteins")
  points(t,solb$p,pch=1,col="red")
  points(sol[,1],sol[,2],col="blue",pch=4)
  par(new=TRUE)
  theo_exact<-ode(y=p0,times=sol[,1],func=fmu,parms=c(ks=el[["ks"]],kd=el[["kd"]]))
  eq2 = paste0("ks=",round(el[["ks"]],3), " kd=",round(el[["kd"]],3)," Fit: ks = ", round(coef(fit)[1],4)," kd=",round(coef(fit)[2],4))
  plot(theo_exact[,1],theo_exact[,2],type="l",xlab="t (days)", ylab="Proteins",ylim=c(0,14),sub=eq2,col="blue")
  legend("top",c("fitted curve","noisified points","exact points","exact curve"),lty=c(1,NA,NA,1),pch=c(NA,1,1,NA),col =c("green","red","blue","blue"),cex=0.5)
  el[["theo"]]<-theo
  el[["theo_exact"]]<-theo_exact
  return(el)
  },error=function(e) NULL)
})
dev.off()




Noisify <- function(data,bruit) {
  # bruit uniforme
  noise <- runif(length(data), -bruit, bruit)
  # bruit suivant une loi normale de variance bruit
  # noise <- rnorm(length(data),0,bruit)
  hist(noise)
  noisified <- data+noise
  while (is_empty(noisified[noisified<0])==FALSE){
    # si la valeur bruit?e est n?gative, on diminue le bruit
    noisified[noisified<0]=data[noisified<0]+runif(length(noisified[noisified<0]), -data[noisified<0], data[noisified<0])
    #noisified[noisified<0]=data[noisified<0]+rnorm(length(noisified[noisified<0]), 0, bruit])
  }
  return(noisified)
}
