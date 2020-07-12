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
plot(t,w,main="poids kiwi")
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
plot(t,mu(t,valp),"l",ylim = c(0,100))
lines(t,sigmoide(t,valp),"l")
##############################################################
# choix rna courbe
rna<-function(t) {
  #out=-10/55^2*(t-55/2)^2+2.5
  out=2*exp(-(t-55/2)^2/500)
  return(out)}

# rna<-function(t) {
#   #out=-10/55^2*(t-55/2)^2+2.5
#   out=2*exp(-t/5)-1/70*t+11/14
#   return(out)}
# rna<-function(t) {
#   #out=-10/55^2*(t-55/2)^2+2.5
#   #out=0.003124*t^2-0.217440*t+4.009098
#   out=0.003559*t^2-0.233527*t+4.078663
#   return(out)}
dt=7
t1=seq(0,max(t),12)
ti=sort(c(t1,t1,t1))
test_rna<-rna(ti)
coef=polyfit(ti,log(test_rna),3)
# g?n?ration des donn?es par r?solution ode


plot(ti,test_rna,"l",col="green")
lines(ti,exp(polyval(coef,ti)),"l",col="red")
# cr?ation donn?es avec prise en compte de mu
##############################################################
fmu=function(t,p,par){
  dp=par["ks"]*exp(polyval(coef,t))-(par["kd"]+mu(t,valp))*p
  return(list(dp))}
p0=3
ks_vrai=5e-6
kd_vrai=5e-3
# parInit<-list("start_prot"=p0,"ks"=ks_vrai,"kd"=kd_vrai)
# parInit2<-list("start_prot"=p0,"ks"=ks_vrai,"kd"=kd_vrai*0.1)
# parInit3<-list("start_prot"=p0,"ks"=ks_vrai*5,"kd"=kd_vrai*5)
# parInit4<-list("start_prot"=p0,"ks"=ks_vrai*10,"kd"=kd_vrai*10)
# listPar<-list("par1"=parInit,"par2"=parInit2,"par3"=parInit3,"par4"=parInit4)
listPar2<-list()
colnames(listPar2)<-c("start_prot","ks","kd")
cont<-1
for (i in seq(1,1000,20)){
  newlist<-list("start_prot"=p0,"ks"=ks_vrai*i,"kd"=kd_vrai*i)
  listPar2[[cont]]=newlist
  cont<-cont+1
}
# listPar<-unlist(listPar2, recursive=F)
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
# cl1 <- makeCluster(detectCores() - 1)
# clusterEvalQ(cl1, {
#   ## set up each worker.  Could also use clusterExport()
#   library(deSolve)
#   library(rlang)# utilisation de is empty
#   library(pracma)
#   library(minpack.lm)
#   NULL
# })
# clusterExport(cl1,c("p0","fmu","Noisify","mu"))
# res_list<-lapply(listPar2, function(el){
#     res<-list()
#     res[["ks_kd"]]=paste0(el[["ks"]],el[["kd"]])
#     sol=ode(y=p0,times=ti,func=fmu,parms=c(ks=el[["ks"]],kd=el[["kd"]]))
#     solb<-matrix(nrow=length(sol[,2]),ncol=2)
#     solb[,1]=sol[,1]
#     bruit=min(sqrt(mean(sol[,2])/6),0.5)
#     val_prot<-Noisify(sol[,2],bruit)
#     solb[,2]=val_prot
#     res[["Noisy sol"]]=solb[,2]
#     return(res)
# })
# res_list<-matrix(nrow = length(ti),ncol = length(listPar2))
# names_col<-c()
# for (i in seq(1,length(listPar2))){
#   sol=ode(y=p0,times=ti,func=fmu,parms=c(ks=listPar2[[i]][["ks"]],kd=listPar2[[i]][["kd"]]))
#   solb<-matrix(nrow=length(sol[,2]),ncol=2)
#   solb[,1]=sol[,1]
#   bruit=min(sqrt(mean(sol[,2])/6),0.5)
#   val_prot<-Noisify(sol[,2],bruit)
#   solb[,2]=val_prot
#   res_list[,i]<-solb[,2]
#   names_col<-append(names_col,paste0("Ks_",listPar2[[i]][["ks"]],"Kd_",kd=listPar2[[i]][["kd"]],sep=""))
# }
# colnames(res_list)<-names_col
# SolB<-read.table("../../test_prot_plusieursKsKd.txt",skip = 1)
# SolB<-SolB[,-1]
load("bruitReal.RData")
SolB<-res_list
for (i in seq(2,2)){
  pdf(sprintf("Fake_data_test_%d_expPoly.pdf",i))
  for (j in seq(1,ncol(SolB))){
    el<-listPar2[[j]]
  # res_list<-lapply(listPar, function(el){
    tryCatch({
    res<-list()
    res[["ks"]]=el[["ks"]]
    res[["kd"]]<-el[["kd"]]
    # dt=7
    # t1=seq(0,55,dt)
    # ti=sort(c(t1,t1,t1))
    sol=ode(y=p0,times=ti,func=fmu,parms=c(ks=el[["ks"]],kd=el[["kd"]]),method = "ode45")
    # solb<-matrix(nrow=length(sol[,2]),ncol=2)
    # solb[,1]=sol[,1]
    # bruit=min( sqrt(mean(sol[,2])/6),0.5)
    # val_prot<-Noisify(sol[,2],bruit)
    # solb[,2]=val_prot
    # res[["Noisy sol"]]=solb
    # res[["Original sol"]]=sol
    # resolution pour ti=0,7,14,...
    #sol=ode(y=p0,times=t,func=f,parms=c(ks=0.31,kd=0.17))
    # plot(solb,ylim=c(0,max(sol[,2])),col="green")
    # points(sol,ylim=c(0,max(sol[,2])),col="blue")
    # curve(rna(x),0,55,col="blue",add=TRUE)
    # legend("topleft",c("pts bruit?s","points exacts","rna"),lty=c(NA,NA,1),pch=c(1,1,NA),col =c("green","blue","blue"),cex=0.5)
    solb=data.frame(t=ti,p=SolB[,j])
    names(solb) <- c("t", "p")
    # attach(solb)
    
    
    ##################################################################"
    # nls pour retrouver ks et kd
    ##################################################################
    
    #nls.control(maxiter = 1000, tol = 1e-03, minFactor = 1/1024,
    #printEval = FALSE, warnOnly = T)
    #ksmin=ks_vrai*0.25#
    ksmin=4.5e-3*mean(test_rna)
    ksmax=1440*mean(test_rna)
    #kdmin=kd_vrai*0.25
    kdmin=4.5e-3
    kdmax=24
    lb=c(ksmin,kdmin,0);ub=c(ksmax,kdmax,2*p0)
    # starting value de kd: max(...) to avoid negative value
    kdinit= el[["kd"]]#max(log(solb[1,2]/solb[4,2])/7,0.001)
    # starting value ks
    ksinit=el[["ks"]] #max(solb[4,2]/solb[4,1]*kdinit,0.001)
    initialisation=list(ks=ksinit,kd=kdinit,p00=mean(solb[1:3,2]))
    fit<-nls(solb$p~ode(y=p00,times=solb[,1],func=fmu,parms=c(ks=ks,kd=kd))[,2],data=list(solb),start=initialisation,trace=TRUE,algorithm="port",control=list(minFactor=1e-6,maxiter=1000,warnOnly=T,tol=sqrt(.Machine$double.eps)),lower=lb,upper=ub)
    browser()
    fit2<-nlsLM(solb$p~ode(y=p00,times=solb[,1],func=fmu,parms=c(ks=ks,kd=kd)),data=list(solb),start=initialisation,lower=lb,upper=ub,trace=TRUE)
    # summary(fit)
    # summary(fit2)
    #fit<-nls(solb$p~ode(y=coef(fit)[3],times=sol[,1],func=f,parms=c(ks=ks,kd=kd))[,2],data=list(solb),start=coef(fit),algorithm="port",control=list(minFactor=1e-5,maxiter=1000),lower=lb,upper=ub)
    res[["coefs"]]=coef(fit)
    theo<-ode(y=coef(fit)[3],times=solb[,1],func=fmu,parms=coef(fit))
    hauteur=max(max(sol[,2]),max(solb[,2]))
    plot(theo[,1],theo[,2],type="l",main="Fit Test ", col="green",ylim=c(0,hauteur+1),xlab="t (days)", ylab="Proteins")
    points(solb$t,solb$p,pch=1,col="red")
    points(sol[,1],sol[,2],col="blue",pch=4)
    par(new=TRUE)
    theo_exact<-ode(y=p0,times=sol[,1],func=fmu,parms=c(ks=res[["ks"]],kd=res[["kd"]]))
    eq2 = paste0("ks=",res[["ks"]], " kd=",res[["kd"]]," Fit: ks = ", coef(fit)[1]," kd=",coef(fit)[2])
    plot(theo_exact[,1],theo_exact[,2],type="l",xlab="t (days)", ylab="Proteins",ylim=c(0,hauteur+1),sub=eq2,col="blue")
    legend("top",c("fitted curve","noisified points","exact points","exact curve"),lty=c(1,NA,NA,1),pch=c(NA,1,1,NA),col =c("green","red","blue","blue"),cex=0.5)
    res[["theo"]]<-theo
    res[["theo_exact"]]<-theo_exact
    # return(res)
    },error=function(e){print("Error in nls fitting") 
      return(NULL)})
  }#)
}
dev.off()



# Noisify <- function(data,bruit) {
#   # bruit uniforme
#   noise <- runif(length(data), -bruit, bruit)
#   # bruit suivant une loi normale de variance bruit
#   # noise <- rnorm(length(data),0,bruit)
#   hist(noise)
#   noisified <- data+noise
#   while (is_empty(noisified[noisified<0])==FALSE){
#     # si la valeur bruit?e est n?gative, on diminue le bruit
#     noisified[noisified<0]=data[noisified<0]+runif(length(noisified[noisified<0]), -data[noisified<0], data[noisified<0])
#     #noisified[noisified<0]=data[noisified<0]+rnorm(length(noisified[noisified<0]), 0, bruit])
#   }
#   return(noisified)
# }
