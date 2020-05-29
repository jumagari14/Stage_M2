# install.packages("minpack.lm")
library("deSolve")
library(minpack.lm)
require("readxl")
require("reshape2")
library(dplyr)
library(pracma)
library(ggplot2)
library(readr)
# setwd("D:/Stage M2/Stage_M2/model/")
setwd("/media/juanma/JUANMA/Stage M2/Stage_M2/model/")

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

# source("input.r")
stageDensity_v2<-function(mrna_data,proteo_data,days_pca,lab_y,lab_x){
  l<-list()
  l[[1]]<-combineGraphs(mrna_data,proteo_data,"All",moyenne = T)
  l[[1]]<-l[[1]]+theme(axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), 
                       axis.ticks.x= element_blank(), 
                       plot.margin= unit(c(0, 1, 0.5, 0.5), "lines"))+scale_x_continuous(trans = "log10",limits = c(1e-6,3e7))
  cont<-1
  unique_dpa<-unique(days_pca)
  for (i in seq(1,ncol(mrna_data),3)){
    cont<-cont+1
    l[[cont]]<-combineGraphs(mrna_data[,i:i+2],proteo_data[,seq(i,i+2)],unique_dpa[cont-1],T)
    
    if (i==ncol(mrna_data)-2){
      l[[cont]]<-l[[cont]]+theme(axis.text.x=element_text(colour="black"),axis.ticks.x=element_line(colour="black"),plot.margin= unit(c(0, 1, 0.5, 0.5), "lines"))+scale_x_continuous(trans = "log10",limits = c(1e-6,3e7))
    }
    else{
      l[[cont]]<-l[[cont]]+theme(axis.title.x = element_blank(), 
                                 axis.text.x = element_blank(), 
                                 axis.ticks.x= element_blank(), 
                                 plot.margin= unit(c(0, 1, 0.5, 0.5), "lines"))+scale_x_continuous(trans = "log10",limits = c(1e-6,3e7))
    }
  }
  # prueba1<-lapply(l, ggplotGrob)
  # prueba2<-do.call(rbind,c(prueba1,size="last"))
  # grid.draw(prueba2)
  plot_grid(plotlist = l,ncol = 1,align = "hv",hjust=1,vjust=1)
  ggarrange(plots = l,ncol=1,bottom = textGrob(lab_x,hjust=0,vjust=0),left =textGrob(lab_y,rot=90,vjust=1))
  # grid.arrange(arrangeGrob(grobs = l,
  #                          nnow = 3,
  #                          left=textGrob(lab_y,rot=90,vjust=1),
  #                          bottom = textGrob(lab_x,hjust=0,vjust=0)),
  #              nrow=1)
}
combineGraphs<-function(mrna_data,proteo_data,annotation,moyenne){
  g<-plotDensity(t(mrna_data),"","","",moyenne)
  if (moyenne){
    melted_pro<-melt(rowMeans(proteo_data,na.rm = T))
  }
  else{
    melted_pro<-melt(proteo_data)
  }
  med_pro<-median(melted_pro$value[melted_pro$value>0],na.rm = T)
  g1<-g+geom_histogram(data = melted_pro,aes(value),bins=30)+geom_vline(xintercept=med_pro,linetype="dotdash",color="darkblue")
  x_max<-10^(ggplot_build(g1)$layout$panel_scales_x[[1]]$range$range)[2]
  x_min<-10^(ggplot_build(g1)$layout$panel_scales_x[[1]]$range$range)[1]
  coord_x_1<-x_max-(x_max-med_pro)*0.3
  coord_x_2<-x_min*0.2+x_min
  y_max<-(ggplot_build(g1)$layout$panel_scales_y[[1]]$range$range)[2]
  g1+annotate("text",x=coord_x_1,y=y_max,hjust=1,vjust=1,label=format(round(med_pro,2),nsmall = 2),color="darkblue")+annotate("text",x=coord_x_2,y=y_max,hjust=0,vjust=1,label=annotation,fontface=2)
}
fitPoids<-function(t,poids,method,dpa_analyse){
  verhulst<-y~(par2*par3)/(par3+(par2-par3)*exp(-par1*t))
  dverhulst_form<-y~r*y*(1-y/K)
  gompertz<-y~par2*exp(log(par3/par2)*exp(-par1*t)) ## ???
  logisgen<-y~par1*(1+exp((par2-t)/par3)^(-1/par4))
  empirique<-y~exp(par1*(t-par3)/(par2+t-par3))
  simpl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))
  doubl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))+par5/(1+exp(-par6*(t-par7)))
  g<-ggplot()+geom_point(aes(x=unlist(t),y=unlist(poids)))
  if (method=="verhulst"){
    ver_fit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=0.1,"par2"=100,"par3"=1),formula_fit  = verhulst,ub=NULL,lb=NULL)
    g<-g+geom_line(aes(x=unlist(t),y=fitted(ver_fit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(ver_fit), names(coef(ver_fit)))
    mu.list <- lapply(mu.list, unname)
    final.form<-ver_fit
  }
  if (method=="gompertz"){
    gom_fit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=0.065,"par2"=114.39,"par3"=0.52),formula_fit = gompertz,ub=NULL,lb=NULL)
    g<-g+geom_line(aes(x=unlist(t),y=fitted(gom_fit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(gom_fit), names(coef(gom_fit)))
    final.form<-gom_fit
  }
  if (method=="empirique"){
    emp_fit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=5.38,"par2"=8,"par3"=7),formula_fit = empirique,ub=NULL,lb=NULL)
    g<-g+geom_line(aes(x=unlist(t),y=fitted(emp_fit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(emp_fit), names(coef(emp_fit)))
    mu.list <- lapply(mu.list, unname)
    final.form<-emp_fit
  }
  if (method=="contois"){
    parList<-c(r=5,K=2,R=100)
    y0<-1
    dev_poids<-odeFitting(as.vector(t),y0,parlist = parList,formula = contois)
    print(dev_poids)
    # cont_fit<-linFitting(as.vector(t),y0,parlist = parList,formula_fit = odeFitting,ub=c(Inf,Inf,130,2),lb=c(0,0,90,0.05))
    # print(fitted(cont_fit))
  }
  if (method=="double_sig"){
    db_sigFit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=40,"par2"=1,"par3"=20,"par4"=1,"par5"=100,"par6"=2,"par7"=45),formula_fit = doubl_sig,ub=c(70,Inf,25,2,110,Inf,55),lb=c(20,0,15,0.1,60,0,40))
    g<-g+geom_line(aes(x=unlist(t),y=fitted(db_sigFit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(db_sigFit), names(coef(db_sigFit)))
    mu.list <- lapply(mu.list, unname)
    final.form<-db_sigFit
    # plot(t,poids)
    # lines(t,fitted(db_sigFit),col=2,lwd=2)
    
  }
  if (method=="log_poly"){
    # poly_model<-lm(log(unlist(poids))~poly(unlist(t),3,raw = T))
    # new_y<-predict.lm(poly_model,data.frame(t))
    wp3<-polyfit(unlist(t),log(unlist(poids)),3)
    new_y<-polyval(wp3,unlist(t))
    new_y<-exp(new_y)
    g<-g+geom_line(aes(x=unlist(t),y=new_y))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- wp3
  }
  print(g)
  return(mu(t,method,mu.list,final.form,dpa_analyse))
# 
#   data_prueba<-select_if(total_data[1:5,],is.numeric)
#   t<-t(as.matrix(seq(1,27)))
#   y<-as.matrix(data_prueba[2,1:27])
#   
#   fitting<-nlsLM(formula =formula_fit,start = parlist,data = data.frame("t"=t,"y"=y))
} 

linFitting<-function(t,y,parlist,formula_fit,ub,lb){
  # colnames(t)<-"t"
  # colnames(y)<-"y"
  fitting<-nlsLM(formula =formula_fit,start = parlist,data = as.data.frame(list("t"=t,"y"=y),col.names = c("t","y")),upper = ub,lower = lb)
  return(fitting)
}
odeFitting<-function(t,y,parlist,formula){
  fitting<-ode(y=c(y=y),t,func = formula,parms = parlist,method = "ode45")
  return(fitting)
}
# plot(t,y)
# lines(t,fitted(fitting),col=2,lwd=2)
loadData<-function(data,trans_sheet,prot_sheet,poids){
  if (grepl("xls[x]*",data,perl = T)){
    data_trans<-read_xlsx(path = data,sheet = trans_sheet)
    data_prot<-read_xlsx(path = data,sheet =prot_sheet)
    colnames(data_prot)[1]<-"Protein"
    colnames(data_trans)[1]<-"Transcrit"
    t<-round(as.numeric(colnames(data_prot)[-which(is.na(as.numeric(as.character(colnames(data_prot)))))]))
    
    total_data<-merge(data_trans,data_prot)
    
    lista<-vector("list",nrow(data_trans))
    for (i in seq(1,nrow(total_data))){
      lista[[i]]<-list("Protein_ID"=total_data[i,"Protein"],"Transcrit_ID"=total_data[i,"Transcrit"],"Transcrit_val"=as.matrix(total_data[i,3:29]),"Protein_val"=as.matrix(total_data[i,30:ncol(total_data)]),"DPA"=t)
      
    }
    return(list("prot"=data_prot,"mrna"=data_trans,"parse"=lista))
  }
  
  if (grepl(".*csv",data,perl = T)){
    total_data<-read_csv(file = data,col_names = F)
    if (poids) colnames(total_data)<-c("DPA","Poids")
    else colnames(total_data)<-total_data[1,]
    return(total_data)
  }
  else{
    total_data<-read_csv(file = data,col_names = F)
    if (poids) colnames(total_data)<-c("DPA","Poids")
    else colnames(total_data)<-total_data[1,]
    return(total_data)
  }
  
  
  
  
}
contois<-function(time,state,par) {
  r<-par["r"]
  K<-par["K"]
  R<-par["R"]
  y<-state["y"]
  dy<-as.numeric(r)*as.numeric(y)*(1-as.numeric(y)/as.numeric(K))/(as.numeric(K)+(as.numeric(R)-1)*y)
  list(dy)
}
normaMean<-function(proteo_data,mrna_data,ks){
  RNA<-mrna_data/mean(mrna_data,na.rm = T)
  PROT<-proteo_data/mean(proteo_data,na.rm = T)
  ks_norm<-ks*mean(mrna_data,na.rm = T)/mean(proteo_data,na.rm = T)
  return(list("mrna"=RNA,"prot"=PROT,"ks"=ks_norm))
}
mu<-function(dpa,method,parlist,formula_fitting,dpa_analyse){
  if (exists("dpa_analyse")){
    dpa<-dpa_analyse
  }
  if (method=="verhulst"){
    coefs<-coef(formula_fitting)
    y_fit<-fitted(formula_fitting)
    mu_val<-function (par1,par2,y) return(par1*(1-y_fit/par2)) ## t ou y
    val<-mu_val(parlist$par1,parlist$par2,dpa)
    
    
  }
  if (method=="gompertz"){
    coefs<-coef(formula_fitting)
    y_fit<-fitted(formula_fitting)
    mu_val<-function (par1,par2,y) return(par1*log(par2/y_fit))
    val<-mu_val(parlist$par1,parlist$par2,dpa)
  }
  if (method=="contois"){
    y0<-parlist$y
    dev_poids<-odeFitting(as.vector(t),y0,parlist = parList,formula = contois)
    y_val<-dev_poids[,2]
    val<-parlist$par1*(1-y_val/parlist$par2)/(parlist$par2+(parlist$par3-1)*y_val)
  }
  if (method=="empirique"){
    val<-parlist$par1*parlist$par2/(parlist$par2+t-parlist$par3)^2
  }
  if (method=="log_poly"){
    d_par<-polyder(parlist)
    val<-polyval(d_par,dpa)/std(dpa)  
    browser()
    ## Acabar...
  }
  if (method=="double_sig"){
    val<-parlist$par2*exp(-parlist$par2*(parlist$par2*(dpa-parlist$par3)))/(1+exp(-parlist$par2*(dpa-parlist$par3)))+parlist$par6*exp(-parlist$par6*(dpa-parlist$par7))/(1+exp(-parlist$par2*(dpa-parlist$par3)))
  }
  return(data.frame("t"=unique(dpa),"Taux"=unique(val)))
}
fit_testRNA<-function(dpa,mrna,fitR){
  model3<-polyfit(dpa,mrna,3)
  model6<-polyfit(dpa,mrna,6)
  model_log<-polyfit(dpa,log(mrna),3)
  
}

# test_data<-lista[[30]]
# 
poids<-read_csv("poids_test.csv",col_names=c("DPA","Poids"))
poids_kiwi<-read_xlsx("Kiwi_FW.xlsx",sheet = "Kiwifruit")
per_dpa<-days_kiwi<-rep(c(0,13,26,39,55,76,118,179,222), each = 3)
test_data<-loadData(data = "Paires_mrna_prot_kiwi_nouvMW.xlsx",trans_sheet = "Transcrits",prot_sheet = "Proteines",F)
test_list<-test_data$parse
taux_croi<-fitPoids(poids_kiwi$DPA,poids_kiwi$Weight_g,"double_sig",per_dpa)
