install.packages("minpack.lm")
library("deSolve")
library(minpack.lm)
require("readxl")
require("reshape2")
library(dplyr)
setwd("/media/juanma/JUANMA/Stage M2/Stage_M2/model/")

verhulst<-nlsLM(formula=y~(par2*par3)/(par3+(par2-par3)*exp(-par1*t)),start=list("par1"=0.1,"par2"=100,"par3"=1),data = data.frame("t"=seq(0,100),"y"=rnorm(101,3,1)))
verhulst<-y~(par2*par3)/(par3+(par2-par3)*exp(-par1*t))
dverhulst_form<-y~r*y*(1-y/K)
gompertz<-y~par2*exp(log(par3/par2)*exp(-par1*t)) ## ???
contois<-y~r*y1*(1-y1/K)/(L+(R-1)*y1)
logisgen<-y~par1*(1+exp((par2-t)/par3)^(-1/par4))
empirique<-y~exp(par1*(t-par3)/(par2+t-par3))

simpl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))
doubl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))+par5/(1+exp(-par6*(t-par7)))

data_prueba<-select_if(total_data[1:5,],is.numeric)
t<-t(as.matrix(seq(1,27)))
y<-as.matrix(data_prueba[2,1:27])

fitting<-nlsLM(formula =doubl_sig,start = list("par1"=110,"par2"=0.1,"par3"=30,"par4"=1),data = data.frame("t"=t,"y"=y))
# plot(t,y)
# lines(t,fitted(fitting),col=2,lwd=2)