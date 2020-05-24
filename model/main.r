# install.packages("minpack.lm")
library("deSolve")
library(minpack.lm)
require("readxl")
require("reshape2")
library(dplyr)
setwd("/media/juanma/JUANMA/Stage M2/Stage_M2/model/")

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
  if (method=="contois"){
    
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

# test_data<-lista[[30]]
# 
# fitPoids(test_data$DPA,test_data$Transcrit_val,"empirique")
