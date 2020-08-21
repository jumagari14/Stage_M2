library("readxl")
library("xlsx")
library("gplots")
library("reshape2")
library("vegan")
library("ggplot2")
library(dplyr)
library("grid")
library("goeveg")
library(gridExtra)
library(cowplot)
library(egg)
setwd("D:/Stage M2/Analyse_stats")
# setwd("/media/juanma/JUANMA/Stage M2/Analyse_stats/")
source("functions.r")

days_kiwi<-rep(c("0 DPA", "13 DPA", "26 DPA", "39 DPA", "55 DPA", "76 DPA", "118 DPA", "179 DPA", "222 DPA"), each = 3)

kiwi_trans<-read.csv("Conc_Red5_filt_corrigé.csv",check.names = F)
colnames(kiwi_trans)[1]<-"Transcrit"

kiwi_prot<-read.csv("Proteo/Prot_conc_modelAc_sansUPS_nouvMW.csv",check.names = F)
kiwi_prot$protein<-NULL
colnames(kiwi_prot)[1]<-"Protein"
kiwi_prot<-kiwi_prot[,c(1,2,3,4,8,9,10,17,18,19,20,21,22,23,24,25,26,27,28,5,6,7,11,12,13,14,15,16)]
kiwi_prot<-kiwi_prot[apply(select_if(kiwi_prot,is.numeric),1,function(x) !all(is.na(x))),]
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


trans_prot<-read.csv("gene_prot_mrna_red5.csv")
colnames(trans_prot)

merge1<-merge(kiwi_trans,trans_prot,by="Transcrit")
merge2<-merge(merge1,kiwi_prot,by="Protein")
merge2[,c("Gene","Accession")]<-NULL

prot_analyse<-merge2[,-c(3:29)]
rownames(prot_analyse)<-prot_analyse$Protein
prot_analyse2<-prot_analyse
prot_analyse$Protein<-NULL
trans_analyse<-merge2[,1:29]
rownames(trans_analyse)<-trans_analyse$Transcrit
trans_analyse2<-trans_analyse
trans_analyse$Transcrit<-NULL

colnames(prot_analyse2)[c(3:29)]<-days_kiwi
colnames(trans_analyse2)[c(3:29)]<-days_kiwi

colnames(prot_analyse)[c(2:28)]<-days_kiwi
colnames(trans_analyse)[c(2:28)]<-days_kiwi

prot_dat<-melt(prot_analyse[15,-1])
trans_dat<-melt(trans_analyse[15,-1])
ligne<-merge(prot_dat,trans_dat,by="variable")

options(repr.plot.width = 25, repr.plot.height = 20)
# pdf("correlation_scatterplot.pdf",width = 25,height=20)
lista<-list()
cont<-1
melt1<-melt(prot_analyse2)
melt2<-melt(trans_analyse2)
melt_def<-merge(melt1,melt2,by = "Transcrit")
melt_def<-melt_def[,c("Transcrit","variable.x","value.x","value.y")]
agregado<-aggregate(melt_def[,3:4],by=list(melt_def$Transcrit),FUN = function(x) mean(x,na.rm = T))
melt_def$value.x<-log10(melt_def$value.x)
melt_def$value.y<-log10(melt_def$value.y)
agregado$value.x<-log10(agregado$value.x)
agregado$value.y<-log10(agregado$value.y)
melt_def[which(is.nan(melt_def$value.x)),]<-NA
melt_def[which(is.infinite(melt_def$value.y)),]<-NA
agregado[which(is.nan(agregado$value.x)),]<-NA
agregado[which(is.infinite(agregado$value.y)),]<-NA
corre<-cor(agregado$value.x,agregado$value.y,use = "complete.obs")
g<-ggplot(data = agregado,aes(x=value.y,y=value.x))+geom_point()+theme
x_max<-(ggplot_build(g)$layout$panel_scales_x[[1]]$range$range)[2]
coord_x<-x_max*0.8
y_min<-(ggplot_build(g)$layout$panel_scales_y[[1]]$range$range)[1]
lista[[1]]<-g+annotate("text",x=coord_x,y=y_min,hjust=1,vjust=1,label=format(round(corre,3),nsmall = 3),color="red")+geom_vline(xintercept = 0,linetype="dashed")+geom_hline(yintercept = 0,linetype="dashed")+xlab("Transcript (fmol per gFW)")+ylab("Protein (fmol per gFW)")
cont<-cont+1
for (i in (seq(3,ncol(prot_analyse2),3))){
  extr_prot<-prot_analyse2[,c(1,2,seq(i,i+2))]
  extr_trans<-trans_analyse2[,c(1,2,seq(i,i+2))]
  melt1<-melt(extr_prot)
  melt2<-melt(extr_trans)
  melt_def<-merge(melt1,melt2,by="Transcrit")
  melt_def<-melt_def[,c("Transcrit","variable.x","value.x","value.y")]
  agregado<-aggregate(melt_def[,3:4],by=list(melt_def$Transcrit,melt_def$variable.x),FUN = function(x) mean(x,na.rm = T))
  melt_def$value.x<-log10(melt_def$value.x)
  melt_def$value.y<-log10(melt_def$value.y)
  agregado$value.x<-log10(agregado$value.x)
  agregado$value.y<-log10(agregado$value.y)
  melt_def[which(is.nan(melt_def$value.x)),]<-NA
  melt_def[which(is.infinite(melt_def$value.y)),]<-NA
  agregado[which(is.nan(agregado$value.x)),]<-NA
  agregado[which(is.infinite(agregado$value.y)),]<-NA
  corre<-cor(agregado$value.x,agregado$value.y,use = "complete.obs")
  g<-ggplot(data = agregado,aes(x=value.y,y=value.x))+geom_point()+theme+xlab("")+scale_y_continuous(limits = c(-1,10))+ylab("")
  x_max<-(ggplot_build(g)$layout$panel_scales_x[[1]]$range$range)[2]
  coord_x<-x_max*0.8
  y_min<-(ggplot_build(g)$layout$panel_scales_y[[1]]$range$range)[1]
  lista[[cont]]<-g+annotate("text",x=coord_x,y=y_min,hjust=1,vjust=1,label=format(round(corre,3),nsmall = 3),color="red")+geom_vline(xintercept = 0,linetype="dashed")+geom_hline(yintercept = 0,linetype="dashed")
  cont<-cont+1
}
m1<-matrix(rep(1,9),nrow = 1)
m2<-matrix(seq(2,10),nrow = 1)
m<-rbind(m1,m2)
# plot_grid(plotlist = lista,ncol = 9,align = "hv",hjust=1,vjust=1)
# ggarrange(plots = lista,ncol=9,bottom = textGrob("",hjust=0,vjust=0),left =textGrob("",rot=90,vjust=1))
grid.arrange(grobs=lista,nrow=2,ncol=9,layout_matrix=m)
# dev.off()

kiwi_prot[-which(kiwi_prot$Protein %in% rownames(prot_analyse)),]->manque
manque_trans<-trans_prot[trans_prot$Protein %in% manque$Protein,c("Transcrit","Protein")]


# 
# write.xlsx(prot_analyse, file="Paires_mrna_prot_kiwi_nouvMW.xlsx", sheetName="Proteines", row.names=T)
# write.xlsx(trans_analyse, file="Paires_mrna_prot_kiwi_nouvMW.xlsx", sheetName="Transcrits",append=T, row.names=T)

## Compare Density 

rownames(kiwi_trans)<-kiwi_trans$Transcrit
kiwi_trans$Transcrit<-NULL
trans_analyse_clean<-select_if(trans_analyse,is.numeric)

# print(plotDensity_v2(kiwi_trans,trans_analyse_clean,""))
# pdf("Analyse_paires_nouvelleMW.pdf",width = 17,height = 22)
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
plotDensity_v2(kiwi_trans,trans_analyse_clean,"","Concentartion (fmol per gFW)","Number of variables")
stageDensity_v2(trans_analyse_clean,select_if(prot_analyse,is.numeric),days_kiwi,"Number of variables","Concentration (fmol per gFW)")
# dev.off()



