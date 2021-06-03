list.of.packages <- c("deSolve", "minpack.lm","readxl","gdata","reshape2","pracma","ggplot2","dplyr","parallel","readr","rlist","data.table","reader","pbapply","car","ellipse","egg","grid")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos="https://pbil.univ-lyon1.fr/CRAN/")
lapply(list.of.packages,require,character.only=TRUE)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

Noisify <- function(data,bruit) {
  noise<-rnorm(length(data),0,bruit)
  noisified <- data+noise
  noisified[noisified<0]=data[noisified<0]*0
  return(noisified)
}
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
  plot_grid(plotlist = l,ncol = 1,align = "hv",hjust=1,vjust=1)
  ggarrange(plots = l,ncol=1,bottom = textGrob(lab_x,hjust=0,vjust=0),left =textGrob(lab_y,rot=90,vjust=1))

}
combineGraphs<-function(mrna_data,proteo_data,annotation,moyenne){
  g<-plotDensity(t(mrna_data),"","","",moyenne)
  if (moyenne){
    melted_pro<-reshape2::melt(rowMeans(proteo_data,na.rm = T))
  }
  else{
    melted_pro<-reshape2::melt(proteo_data)
  }
  med_pro<-median(melted_pro$value[melted_pro$value>0],na.rm = T)
  g1<-g+geom_histogram(data = melted_pro,aes(value),bins=30)+geom_vline(xintercept=med_pro,linetype="dotdash",color="darkblue")
  x_max<-10^(ggplot_build(g1)$layout$panel_scales_x[[1]]$range$range)[2]
  x_min<-10^(ggplot_build(g1)$layout$panel_scales_x[[1]]$range$range)[1]
  coord_x_1<-x_max-(x_max-med_pro)*0.3
  coord_x_2<-x_min*0.2+x_min
  y_max<-(ggplot_build(g1)$layout$panel_scales_y[[1]]$range$range)[2]
  g1<-g1+annotate("text",x=coord_x_1,y=y_max,hjust=1,vjust=1,label=format(round(med_pro,2),nsmall = 2),color="darkblue")+annotate("text",x=coord_x_2,y=y_max,hjust=0,vjust=1,label=annotation,fontface=2)+xlab(bquote("Concentration "~(fmolgFW^-1)))+ylab("Number of occurences")
}
plotDensity<-function(dataframe, title,lab_x,lab_y,moyenne){
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
  if (moyenne){
    melted<-reshape2::melt(colMeans(dataframe,na.rm = T))
  }
  else{
    melted<-reshape2::melt(dataframe)
  }
  med<-median(melted$value[melted$value>0],na.rm = T)
  g<-ggplot()+geom_histogram(data=melted,aes(value),bins = 30)+scale_x_log10()+xlab(lab_x)+
    ylab(lab_y) +theme+ggtitle(title)+
    geom_vline(xintercept=med,linetype="dotdash",color="red")
  x_min<-10^(ggplot_build(g)$layout$panel_scales_x[[1]]$range$range)[1]
  coord_x<-(med-x_min)*0.3+x_min
  y_max<-(ggplot_build(g)$layout$panel_scales_y[[1]]$range$range)[2]
  r<-g+annotate("text",x=coord_x,y=y_max,hjust=1,vjust=1,label=format(round(med,3),nsmall = 3),color="red")
  return(r)
}
resultsKsKdUI<-function(id){
  ns<-NS(id)
  tagList(htmlOutput(ns("scores")),
          plotOutput(ns("fits")),
          plotOutput(ns("ellip")))
}
resultsKsKd<-function(input,output,session,resvalid,trans_id){
  ns<-session$ns
  pair_ret<<-Filter(function(x) identical(x[["Transcrit_ID"]],trans_id),resvalid)
  output$scores<-renderUI({
    ks1<-paste("ks value",pair_ret[[1]][["SOL"]][["solK"]][2,1]*mean(pair_ret[[1]][["Protein_val"]],na.rm=T)/mean(pair_ret[[1]][["Transcrit_val"]],na.rm = T))
    ks2<-paste("Ks value",pair_ret[[1]][["SOL"]][["solK"]][2,1])
    ks3<-paste("Kd value",pair_ret[[1]][["SOL"]][["solK"]][3,1])
    mes0<-paste("Protein fitting error: ",pair_ret[[1]][["SOL"]][["error"]][["errg"]][1])
    mes1<-paste("Error message: ",pair_ret[[1]][["SOL"]][["error"]][["message"]])
    mes2<-paste("Score of protein fitting error: ",pair_ret[[1]][["SOL"]][["error"]][["score"]])
    mes3<-paste("Optimization message: ",pair_ret[[1]][["SOL"]][["opt_eval"]][["message"]])
    mes4<-paste("Optimization score: ",pair_ret[[1]][["SOL"]][["opt_eval"]][["score"]])
    mes5<-paste("mRNA fitting error: ",pair_ret[[1]][["errorMrna"]])
    mes6<-paste("Weight fitting error: ",pair_ret[[1]][["errorWeight"]])
    HTML(paste(ks2,ks1,ks3,mes0,mes1,mes2,mes3,mes4,mes5,mes6, sep = '<br/>'))
  })
  output$fits<-renderPlot({
    grid.arrange(pair_ret[[1]][["plot_mrna"]],pair_ret[[1]][["SOL"]][["plot_fit_prot"]],ncol=2)
  })
  output$ellip<-renderPlot({
    pair_ret[[1]][["SOL"]][["confEllipsePlot"]]+geom_point(aes(x=pair_ret[[1]][["SOL"]][["solK"]][2,1],y=pair_ret[[1]][["SOL"]][["solK"]][3,1],shape="circle filled"),fill="red")+scale_shape(guide=FALSE)
  })
  
}
paramListInput<-function(id){
  ns<-NS(id)
  tagList(uiOutput(ns("listpar")))
}
paramList<-function(input,output,session,method){
  ns <- session$ns
  if (method=="verhulst"){
    output$listpar<-renderUI({
      tabsetPanel(type = "pills",
                  tabPanel("Parameters",
                           tagList(textInput(ns("par1_sig"),"Enter value of a",value =0.1),
                                   textInput(ns("par2_sig"),"Enter value of b",value = 100),
                                   textInput(ns("par3_sig"),"Enter value of c",value = 1))),
                  tabPanel("Boundaries",
                           tagList(sliderInput(ns("bound_a"), "Upper and lower bounds of a",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_b"), "Upper and lower bounds of b",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_c"), "Upper and lower bounds of c",min = -5000, max = 5000,value = c(0,5000)))))
    })
  }
  if (method=="double_sig"){
    output$listpar<-renderUI({
      tabsetPanel(type="pills",
                  tabPanel("Parameters",
                           tagList(textInput(ns("par1_sig"),"Enter value of a",value = 48),
                                   textInput(ns("par2_sig"),"Enter value of b",value = 0.144),
                                   textInput(ns("par3_sig"),"Enter value of c",value = 35),
                                   textInput(ns("par4_sig"),"Enter value of d",value = 0.4),
                                   textInput(ns("par5_sig"),"Enter value of e",value = 48),
                                   textInput(ns("par6_sig"),"Enter value of f",value = 0.042),
                                   textInput(ns("par7_sig"),"Enter value of g",value = 90))),
                  tabPanel("Boundaries",
                           tagList(sliderInput(ns("bound_a"), "Upper and lower bounds of a",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_b"), "Upper and lower bounds of b",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_c"), "Upper and lower bounds of c",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_d"), "Upper and lower bounds of d",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_e"), "Upper and lower bounds of e",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_f"), "Upper and lower bounds of f",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_g"), "Upper and lower bounds of g",min = -5000, max = 5000,value = c(0,5000)))))
    })
  }
  if (method=="empirique"){
    output$listpar<-renderUI({
      tabsetPanel(type = "pills",
                  tabPanel("Parameters",
                           tagList(textInput(ns("par1_sig"),"Enter value of a",value =5.38),
                                   textInput(ns("par2_sig"),"Enter value of b",value = 8),
                                   textInput(ns("par3_sig"),"Enter value of c",value = 7))),
                  tabPanel("Boundaries",
                           tagList(sliderInput(ns("bound_a"), "Upper and lower bounds of a",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_b"), "Upper and lower bounds of b",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_c"), "Upper and lower bounds of c",min = -5000, max = 5000,value = c(0,5000)))))
      
    })
  }
  if (method=="gompertz"){
    output$listpar<-renderUI({
      tabsetPanel(type="pills",
                  tabPanel("Parameters",
                           tagList(textInput(ns("par1_sig"),"Enter value of a",value =0.065),
                                   textInput(ns("par2_sig"),"Enter value of b",value = 114.39),
                                   textInput(ns("par3_sig"),"Enter value of c",value = 0.52))),
                  tabPanel("Boundaries",
                           tagList(sliderInput(ns("bound_a"), "Upper and lower bounds of a",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_b"), "Upper and lower bounds of b",min = -5000, max = 5000,value = c(0,5000)),
                                   sliderInput(ns("bound_c"), "Upper and lower bounds of c",min = -5000, max = 5000,value = c(0,5000)))))
      
    })
  }
  if (method=="log_poly"){
    output$listpar<-renderUI({
      tagList()})
  }
  
}

getParams<-function(input,method){
  x<-reactiveValuesToList(input)
  if (any(method==c("verhulst","empirique","gompertz"))){
    x[["params-par4_sig"]]<-NULL
    x[["params-par5_sig"]]<-NULL
    x[["params-par6_sig"]]<-NULL
    x[["params-par7_sig"]]<-NULL
    x[["params-bound_d"]]<-NULL
    x[["params-bound_e"]]<-NULL
    x[["params-bound_f"]]<-NULL
    x[["params-bound_g"]]<-NULL
  }
  if (method=="log_poly"){
    x[["params-par1_sig"]]<-NULL
    x[["params-par2_sig"]]<-NULL
    x[["params-par3_sig"]]<-NULL
    x[["params-par4_sig"]]<-NULL
    x[["params-par5_sig"]]<-NULL
    x[["params-par6_sig"]]<-NULL
    x[["params-par7_sig"]]<-NULL
    x[["params-bound_a"]]<-NULL
    x[["prams-bound_b"]]<-NULL
    x[["params-bound_c"]]<-NULL
    x[["params-bound_d"]]<-NULL
    x[["prams-bound_e"]]<-NULL
    x[["params-bound_f"]]<-NULL
    x[["params-bound_g"]]<-NULL
  }
  x_ind<-grep("par[1-9]+",names(x),perl = T)
  newlist<-vector("list",length(x_ind))
  names(newlist)<-names(x[x_ind])
  for (el in names(newlist)){
    newlist[[el]]<-as.numeric(as.character(input[[el]]))
  }
  names(newlist)<-gsub("_sig","",names(newlist))
  names(newlist)<-gsub("params-","",names(newlist))
  
  x_ind<-grep("bound_[a-z]+",names(x),perl = T)
  bounds<-list()
  bounds[["ub"]]<-vector("double",length = length(x_ind))
  bounds[["lb"]]<-vector("double",length = length(x_ind))
  cont<-1
  for (el in names(x[x_ind])){
    bounds[["ub"]][cont]<-as.numeric(as.character(input[[el]][2]))
    bounds[["lb"]][cont]<-as.numeric(as.character(input[[el]][1]))
    cont<-cont+1
  }
  return(list("para"=newlist,"bounds"=bounds))
}

fitPoids<-function(t,poids,method){
  verhulst<-y~(par2*par3)/(par3+(par2-par3)*exp(-par1*t))
  dverhulst_form<-y~r*y*(1-y/K)
  gompertz<-y~par2*exp(log(par3/par2)*exp(-par1*t)) 
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
  if (method=="double_sig"){
    parlist = list("par1"=48,"par2"=0.144,"par3"=35,"par4"=0.4,"par5"=48,"par6"=0.042,"par7"=90)
    db_sigFit<-linFitting(as.vector(t),as.vector(poids),parlist = parlist,formula_fit = doubl_sig,ub=NULL,lb=c(0,0,0,0,0,0,0))
    g<-g+geom_line(aes(x=unlist(t),y=fitted(db_sigFit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(db_sigFit), names(coef(db_sigFit)))
    mu.list <- lapply(mu.list, unname)
    final.form<-db_sigFit
    
  }
  if (method=="log_poly"){
    wp3<-polyfit(unlist(t),log(unlist(poids)),3)
    new_y<-polyval(wp3,unlist(t))
    new_y<-exp(new_y)
    g<-g+geom_line(aes(x=unlist(t),y=new_y))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- wp3
  }
  err<-norm(poids-fitted(final.form),"2")/norm(poids,"2")
  return(list("coefs"=mu.list,"formula"=final.form,"error"=err,"graph"=g))
} 
fitPoids_v2<-function(t,poids,method,listpar,bounds){
  t<-as.numeric(pull(t))
  poids<-as.numeric(pull(poids))
  poids<-poids[!is.na(poids)]
  t<-t[!is.na(t)]
  
  verhulst<-y~(par2*par3)/(par3+(par2-par3)*exp(-par1*t))
  dverhulst_form<-y~r*y*(1-y/K)
  gompertz<-y~par2*exp(log(par3/par2)*exp(-par1*t)) 
  logisgen<-y~par1*(1+exp((par2-t)/par3)^(-1/par4))
  empirique<-y~exp(par1*(t-par3)/(par2+t-par3))
  simpl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))
  doubl_sig<-y~par4+par1/(1+exp(-par2*(t-par3)))+par5/(1+exp(-par6*(t-par7)))
  g<-ggplot()+geom_point(aes(x=unlist(t),y=unlist(poids)),size=3)
  
  
  if (method=="verhulst"){
    ver_fit<-linFitting(as.vector(t),as.vector(poids),parlist = listpar,formula_fit  = verhulst,ub=bounds$ub,lb=bounds$lb)
    g<-g+geom_line(aes(x=unlist(t),y=fitted(ver_fit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(ver_fit), names(coef(ver_fit)))
    mu.list <- lapply(mu.list, unname)
    final.form<-ver_fit
    err<-norm(poids-fitted(ver_fit),"2")/norm(poids,"2")
  }
  if (method=="gompertz"){
    gom_fit<-linFitting(as.vector(t),as.vector(poids),parlist = listpar,formula_fit = gompertz,ub=bounds$ub,lb=bounds$lb)
    g<-g+geom_line(aes(x=unlist(t),y=fitted(gom_fit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(gom_fit), names(coef(gom_fit)))
    final.form<-gom_fit
    err<-norm(poids-fitted(gom_fit),"2")/norm(poids,"2")
  }
  if (method=="empirique"){
    emp_fit<-linFitting(as.vector(t),as.vector(poids),parlist = list("par1"=5.38,"par2"=8,"par3"=7),formula_fit = empirique,ub=bounds$ub,lb=bounds$lb)
    g<-g+geom_line(aes(x=unlist(t),y=fitted(emp_fit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(emp_fit), names(coef(emp_fit)))
    mu.list <- lapply(mu.list, unname)
    final.form<-emp_fit
    err<-norm(poids-fitted(emp_fit),"2")/norm(poids,"2")
  }
  if (method=="double_sig"){
    
    db_sigFit<-linFitting(as.vector(t),as.vector(poids),parlist = listpar,formula_fit = doubl_sig,ub=bounds$ub,lb=bounds$lb)
    
    g<-g+geom_line(aes(x=unlist(t),y=fitted(db_sigFit)))+theme+xlab("DPA")+ylab("Weight (gFW)")
    mu.list <- split(coef(db_sigFit), names(coef(db_sigFit)))
    mu.list <- lapply(mu.list, unname)
    final.form<-db_sigFit
    err<-norm(poids-fitted(db_sigFit),"2")/norm(poids,"2")
    
  }
  return(list("coefs"=mu.list,"formula"=final.form,"graph"=g,"error"=err))
} 

linFitting<-function(t,y,parlist,formula_fit,ub,lb){
  names(t)<-"t"
  names(y)<-"y"
  fitting<-nlsLM(formula =formula_fit,start = parlist,data = as.data.frame(list("t"=as.matrix(t),"y"=as.matrix(y))),upper = ub,lower = lb)
  return(fitting)
}

loadData<-function(data,trans_sheet,prot_sheet,poids){
  if (grepl("xls[x]*",data,perl = T)){
    data_trans<-read_xlsx(path = data,sheet = trans_sheet)
    data_prot<-read_xlsx(path = data,sheet =prot_sheet)
  }
  else{
    data_trans<-read_delim(file = data[1],delim=get.delim(data[1]),col_names = T,na=c("","NA","#E/E","NaN"))
    data_prot<-read_delim(file = data[2],delim=get.delim(data[2]),col_names = T,na=c("","NA","#E/E","NaN"))
  }
  data_prot[,3:ncol(data_prot)]<-sapply(data_prot[,3:ncol(data_prot)], as.numeric)
  data_trans[,3:ncol(data_trans)]<-sapply(data_trans[,3:ncol(data_trans)], as.numeric)
  colnames(data_prot)[1]<-"Protein"
  colnames(data_trans)[1]<-"Transcrit"
  colnames(data_prot)<-gsub("\\.\\.\\.[0-9]*","",colnames(data_prot))
  t<-round(as.numeric(colnames(data_prot)[-which(is.na(as.numeric(as.character(colnames(data_prot)))))]))
  
  total_data<-merge(data_trans,data_prot,by=c("Protein","Transcrit"))
  
  lista<-vector("list",nrow(data_trans))
  for (i in seq(1,nrow(total_data))){
    lista[[i]]<-list("Protein_ID"=total_data[i,"Protein"],"Transcrit_ID"=total_data[i,"Transcrit"],"Transcrit_val"=as.matrix(total_data[i,3:ncol(data_trans)]),"Protein_val"=as.matrix(total_data[i,(ncol(data_trans)+1):ncol(total_data)]),"DPA"=t)
    
  }
  return(list("prot"=data_prot,"mrna"=data_trans,"parse"=lista))
  
  
  
}
normaMean<-function(proteo_data,mrna_data,ks){
  RNA<-mrna_data/mean(mrna_data,na.rm = T)
  PROT<-proteo_data/mean(proteo_data,na.rm = T)
  ks_norm<-ks*mean(mrna_data,na.rm = T)/mean(proteo_data,na.rm = T)
  
  return(list("mrna"=RNA,"prot"=PROT,"ks"=ks_norm))
}
mu<-function(dpa,method,parlist,formula_fitting,dpa_analyse){
  if (method=="verhulst"){
    coefs<-coef(formula_fitting)
    y_fit<-fitted(formula_fitting)
    mu_val<-function (par1,par2,y) return(par1*(1-y_fit/par2)) 
    val<-mu_val(parlist$par1,parlist$par2,dpa)
    
    
  }
  if (method=="gompertz"){
    coefs<-coef(formula_fitting)
    y_fit<-fitted(formula_fitting)
    mu_val<-function (par1,par2,y) return(par1*log(par2/y_fit))
    val<-mu_val(parlist$par1,parlist$par2,dpa)
  }
  if (method=="empirique"){
    val<-parlist$par1*parlist$par2/(parlist$par2+dpa-parlist$par3)^2
  }
  if (method=="log_poly"){
    d_par<-polyder(unlist(parlist,use.names = F))
    val<-polyval(d_par,dpa)/std(dpa)  
    
  }
  if (method=="double_sig"){
    
    val<-parlist$par1*parlist$par2*exp(-parlist$par2*(dpa-parlist$par3))/(1+exp(-parlist$par2*(dpa-parlist$par3)))^2+parlist$par5*parlist$par6*exp(-parlist$par6*(dpa-parlist$par7))/(1+exp(-parlist$par6*(dpa-parlist$par7)))^2
  }
  return(val)
}
fit_testRNA<-function(dpa,mrna,fitR){
  model3<-polyfit(dpa,mrna,3)
  model6<-polyfit(dpa,mrna,6)
  model_log<-polyfit(dpa,log(mrna),3)
  if (fitR=="3_deg"){
    ret<-model3
  }
  if (fitR=="6_deg"){
    ret<-model6
  }
  if (fitR=="3_deg_log"){
    ret<-model_log
  }
  val_fit<-solmRNA(dpa,ret,fitR)
  fitErrMRNA<-norm(mrna-val_fit,"2")/norm(mrna,"2")
  
  return(list("coefs"=ret,"error"=fitErrMRNA))
}
solgss_Borne<-function(dpa,prot_conc,ks_min,ksnorm,algo){
  
  init_prot<-init_conc(dpa,prot_conc)
  if (is.nan(init_prot$init)){
    return(NULL)
  }
  else{
    if (any(is.na(prot_conc))){
      index_nan<-which(is.na(prot_conc))
      dpa<-dpa[-index_nan]
      prot_conc<-prot_conc[-index_nan]
    }
    init<-init_prot$init
    parInit<-list("start_prot"=init,"ks"=ks_min*3,"kd"=ks_min*0.3)
    parInit2<-list("start_prot"=init,"ks"=ks_min,"kd"=ks_min*0.1)
    parInit3<-list("start_prot"=init,"ks"=ks_min*0.05,"kd"=ks_min*0.5)
    parInit4<-list("start_prot"=init,"ks"=ks_min*0.01,"kd"=ks_min*0.1)
    
    if (algo=="LM"){
      parMu<-nlsLM(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit,control=nls.lm.control(ftol=1e-6,maxiter=1000,ptol=1e-6),lower = c(0,0,0))
      parMu2<-nlsLM(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit2,control=nls.lm.control(ftol=1e-6,maxiter=1000,ptol=1e-6),lower = c(0,0,0))
      parMu3<-nlsLM(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit3,control=nls.lm.control(ftol=1e-6,maxiter=1000,ptol=1e-6),lower = c(0,0,0))
      parMu4<-nlsLM(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit4,control=nls.lm.control(ftol=1e-6,maxiter=1000,ptol=1e-6),lower = c(0,0,0))
      
    }
    else if (algo=="Port"){
      parMu<-nls(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit,control=list(minFactor=1e-6,maxiter=1000,tol=1e-6),algorithm = "port",lower = c(0,0,0))
      parMu2<-nls(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit2,control=list(minFactor=1e-6,maxiter=1000,tol=1e-6),algorithm = "port",lower = c(0,0,0))
      parMu3<-nls(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit3,control=list(minFactor=1e-6,maxiter=1000,tol=1e-6),algorithm = "port",lower = c(0,0,0))
      parMu4<-nls(prot_conc~ode(y=start_prot,times = dpa,func = eqDifPrinc,parms = c(ks=ks,kd=kd),method = "ode45")[,2],start = parInit4,control=list(minFactor=1e-6,maxiter=1000,tol=1e-6),algorithm = "port",lower = c(0,0,0))
      
    }
    
    sol1<-fitted(parMu)
    sol2<-fitted(parMu2)
    sol3<-fitted(parMu3)
    sol4<-fitted(parMu4)
    sol<-cbind(sol1,sol2,sol3,sol4)
    parmu<-cbind(coef(parMu),coef(parMu2),coef(parMu3),coef(parMu4))
    resnorm<-c(deviance(parMu),deviance(parMu2),deviance(parMu3),deviance(parMu4))
    
    opt_setp<-evalOptim(parmu)
    error<-sqrt(resnorm[1])/norm(prot_conc,"2")
    errg<-sqrt(resnorm)/norm(prot_conc,"2")
    err_mes<-scoreErreur(error)
    err_mes[["errg"]]<-errg
    model_list<-list("model1"=parMu,"model2"=parMu2,"model3"=parMu3,"model4"=parMu4)
    confEllipse<-confidenceEllipse(parMu,which.coef = c("ks","kd"),fill = T,segments = 50,levels=c(0.9,0.75,0.5))
    bound_circle<-confEllipse[["0.9"]]
    dist_center_axis<-t(t(bound_circle)-c(coef(parMu)["ks"],coef(parMu)["kd"]))
    dist2center<-sqrt(rowSums(dist_center_axis^2))
    max_dis<-abs(dist_center_axis[which(dist2center==max(dist2center)),])
    if ((max_dis[1]<max(coef(parMu)[["ks"]],0.5)) & (max_dis[2]<max(coef(parMu)[["kd"]],0.5))){
      validEllipse<-TRUE
    }
    else{
      validEllipse<-FALSE
    }
    confEllipse2<-do.call(gdata::combine,confEllipse)
    if (any(confEllipse2[,c("x","y")]<0)){
      confEllipsePlot<-ggplot()+geom_polygon(confEllipse2,mapping=aes(x=x,y=y,colour=source,fill=source),alpha=0.3)+theme+coord_cartesian(xlim=(c(0,max(confEllipse2$x))),ylim = (c(0,max(confEllipse2$y))))+ylab("kd")+xlab("Ks")
    } 
    else {
      confEllipsePlot<-ggplot()+geom_polygon(confEllipse2,mapping=aes(x=x,y=y,colour=source,fill=source),alpha=0.3)+theme+ylab("kd")+xlab("Ks")
    }
    
    return(list("solK"=parmu,"sumsq"=resnorm,"opt_eval"=opt_setp,"error"=err_mes,"prot_fit"=sol,"modelList"=model_list,"validEllipse"=validEllipse,"confEllipsePlot"=confEllipsePlot))
  }
  
}

evalOptim<-function(parmu){
  parmu1<-as.vector(parmu[,1])
  parmu2<-as.vector(parmu[,2])
  parmu3<-as.vector(parmu[,3])
  parmu4<-as.vector(parmu[,4])
  
  if (all(pmax(parmu1,parmu2,na.rm = T)>1e-4)){
    ecart12<-abs(parmu1-parmu2)/pmax(parmu1,parmu2,na.rm=T)
  }
  else{
    ecart12<-abs(parmu1-parmu2)
  }
  if (all(pmax(parmu2,parmu3,na.rm = T)>1e-4)){
    ecart23<-abs(parmu2-parmu3)/pmax(parmu2,parmu3,na.rm=T)
  }
  else{
    ecart23<-abs(parmu2-parmu3)
  }
  if (all(pmax(parmu3,parmu4,na.rm = T)>1e-4)){
    ecart34<-abs(parmu3-parmu4)/pmax(parmu3,parmu4,na.rm=T)
  }
  else{
    ecart34<-abs(parmu3-parmu4)
  }
  if (all(ecart12<5e-2) & (all(ecart23<5e-2)) & (all(ecart34<5e-2))){
    mess<-"Optimisation converges to 3 indentical values"
    score<-10
  }
  else if (all(ecart12<1e-1) & (all(ecart23<1e-1)) & (all(ecart34<1e-1))){
    mess<-"Optimisation nearly converges to 3 identical values"
    score<-8
  }
  else if (all(ecart12<5e-2)){
    score<-6
    mess<-"Optimisation converges to 2 indentical values"
  }
  else if (all(ecart23<1e-1)){
    score<-4
    mess<-"Optimisation nearly converges to 2 identical values"
  }
  else{
    score<-1
    mess<-"Optimisation converges to different values"
  }
  return(list("score"=score,"message"=mess))
}
scoreErreur<-function(erreur){
  if (erreur<0.1){
    mess<-"Excellent fit"
  }
  else if ((0.1<erreur) & (erreur<0.15)){
    mess<-"Very accurate fit"
  }
  else if ((0.15<erreur) & (erreur<0.2)){
    mess<-"Accurate fit"
  }
  else if ((0.2<erreur) & (erreur<0.3)){
    mess<-"Pretty accurate fit"
  }
  else if ((0.3<erreur) & (erreur<0.4)){
    mess<-"Average fit"
  }
  else{
    mess<-"Inaccurate fit"
  }
  score<-floor((0.5-erreur)*20)
  return(list("score"=score,"message"=mess))
}
solmRNA<-function(dpa,coef_list,fitR){
  val<-polyval(coef_list,dpa)
  if (fitR=="3_deg_log"){
    val<-exp(val)
  }
  return(val)
}
init_conc<-function(dpa,prot_conc){
  dpa_min<-min(dpa)
  ind<-which(dpa==dpa_min)
  
  if (length(ind)==1){
    p0<-prot_conc[ind]
    pMax<-5*p0
    pMin<-0.1*p0
  }
  else{
    p0<-mean(prot_conc[ind],na.rm = T)
    pMin<-min(prot_conc[ind])-std(prot_conc)*2
    pMax<-max(prot_conc[ind])+std(prot_conc)*2
  }
  return(list("init"=p0,"min"=pMin,"max"=pMax))
}
eqDifPrinc<-function(time,state,par){
  ks<-par["ks"]
  kd<-par["kd"]
  val<-unlist(ks)*solmRNA(time,fittedmrna$coefs,fitR)-(unlist(kd)+mu(dpa=c(time),fitWe,poids_coef,formula_poids,dpa_analyse = NULL)/predict(formula_poids,data.frame(t=time)))*state
  
  return(list(val))
}

plotFitmRNA<-function(dpa,exp_data,fit_data,title){
  
  merg_data<-as.data.frame(cbind(dpa,as.vector(exp_data),fit_data))
  g<-ggplot(data = merg_data,aes(x=dpa))+geom_point(aes(y=V2))+theme+xlab("Time (DPA)")+ylab("mRNA concentration normalized by mean (fmol gFW)")+geom_line(aes(y=fit_data))+ggtitle(title)
  return(g)
}
plotFitProt<-function(dpa,exp_data,fit_data,title){
  if(any(is.na(exp_data))){
    exp_data<-exp_data[-which(is.na(exp_data))]
    dpa<-dpa[-which(is.na(exp_data))]
  }
  merg_data<-as.data.frame(cbind(dpa,exp_data,fit_data))
  melt_data<-reshape2::melt(merg_data[,c("dpa","sol1","sol2","sol3","sol4")],id.vars="dpa")
  g<-ggplot(data = merg_data,aes(x=dpa))+geom_point(aes(y=exp_data))+theme+xlab("Time (DPA)")+ylab("Protein concentration normalized by mean (fmol gFW)")+geom_line(data = melt_data,aes(x=dpa,y=value,group=variable))+ggtitle(title)
  return(g)
}
