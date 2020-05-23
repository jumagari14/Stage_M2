# install.packages("minpack.lm")
library("deSolve")
library(minpack.lm)
require("readxl")
require("reshape2")
library(dplyr)
setwd("/media/juanma/JUANMA/Stage M2/Stage_M2/model/")

source("input.r")

fitPoids<-function(t,poids,method){
  verhulst<-y~(par2*par3)/(par3+(par2-par3)*exp(-par1*t))
  dverhulst_form<-y~r*y*(1-y/K)
  gompertz<-y~par2*exp(log(par3/par2)*exp(-par1*t)) ## ???
  contois<-y~r*y1*(1-y1/K)/(L+(R-1)*y1)
  logisgen<-y~par1*(1+exp((par2-t)/par3)^(-1/par4))
  empirique<-y~exp(par1*(t-par3)/(par2+t-par3))
  simpl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))
  doubl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))+par5/(1+exp(-par6*(t-par7)))
  if (method=="verhulst"){
    ver_fit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=0.1,"par2"=100,"par3"=1),formula = verhulst,ub=NULL,lb=NULL)
    print(ver_fit)
  }
  if (method=="gompertz"){
    gom_fit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=0.065,"par2"=114.39,"par3"=0.52),formula = gompertz,ub=NULL,lb=NULL)
    print(gom_fit)
  }
  if (method=="empirique"){
    emp_fit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=5.38,"par2"=8,"par3"=7),formula = empirique,ub=NULL,lb=NULL)
    print(predict(emp_fit))
  }
  
# 
#   data_prueba<-select_if(total_data[1:5,],is.numeric)
#   t<-t(as.matrix(seq(1,27)))
#   y<-as.matrix(data_prueba[2,1:27])
#   
#   fitting<-nlsLM(formula =formula_fit,start = parlist,data = data.frame("t"=t,"y"=y))
} 

linFitting<-function(t,y,parlist,formula_fit,ub,lb){
  fitting<-nlsLM(formula =formula_fit,start = parlist,data = data.frame(t=t,y=y),upper = ub,lower = lb)
  return(fitting)
}
odeFitting<-function(t,y,parlist,formula){
  fitting<-ode(y,1,func = formula,parms = parlist,method = "ode45")
}
# plot(t,y)
# lines(t,fitted(fitting),col=2,lwd=2)


# fitPoids<-function(xi,yi,opt_func){
#   switch (opt_func,
#     verhulst = verhulstFunc(xi,yi)
#   )
# }

test_data<-lista[[30]]

fitPoids(test_data$DPA,test_data$Transcrit_val,"empirique")
