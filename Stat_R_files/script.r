# Load packages -----------------------------------------------------------
setwd("E:/Stage M2/Analyse_stats")
# setwd("/media/juanma/JUANMA/Stage M2/Analyse_stats/")
library("readxl")
library("gplots")
library("reshape2")
library("vegan")
library("ggplot2")
library("grid")
library(cowplot)
library(egg)
library(ggrepel)
source("functions.r")

# Load Data ---------------------------------------------------------------
donnees_kiwi<-getData("Conc_Red5.csv",",")
days_kiwi<-rep(c("0 DPA", "13 DPA", "26 DPA", "39 DPA", "55 DPA", "76 DPA", "118 DPA", "179 DPA", "222 DPA"), each = 3)

kiwi_data<-donnees_kiwi$format2
kiwi_data_2<-donnees_kiwi$format1

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

# Densities ---------------------------------------------------------------
stageDensity(kiwi_data,days_kiwi,"Concentration (fmol per gFW)","Number of variables",moyenne = T)
plotDensity(kiwi_data,"","Concentration (fmol per gFW)","Number of occurences",T)

# Heatmap -----------------------------------------------------------------
getHeatMap(kiwi_data_2,"Correlation heatmap for kiwi Red5")
clusteredHeatmap(kiwi_data_2,"Heatmap of conc for kiwi")
# #---------------------------------------------------------------------

getPCA(kiwi_data,days_kiwi,"ACP of kiwi concentration")
stage_kiwi<-aggregate(kiwi_data,list(rep(1:(nrow(kiwi_data)%/%3+1),each=3,len=nrow(kiwi_data))),function(x) mean(x,na.rm = T))[-1]
rownames(stage_kiwi)<-unique(days_kiwi)
getPCA(stage_kiwi,unique(days_kiwi),"PCA of mRNA concentration")

var_kiwi<-apply(stage_kiwi,1,function(x) cv(x)*100)
data<-data.frame(var_kiwi)
colnames(data)<-c("Kiwi")
data<-melt(data)
data$stade<-as.factor(rep(seq(1:9),3))
colnames(data)<-c("Espèce","value","stade")
ggplot(data, aes(fill=Esp?ce, y=value, x=stade)) +
  geom_bar(position="dodge", stat="identity")+theme+scale_fill_manual(values=c("yellow","green","red"))+ylab("CV (%)")



# Annot per stage ---------------------------------------------------------

annot_kiwi<-read.csv("Prot_red5.csv",header = T,sep = "\t",check.names = F)
annot_kiwi<-annot_kiwi[,-1]

gene_mrna<-read.csv("gene_mrna_Red5.csv",check.names = F)

print(medianAnno_trans(annot_kiwi,kiwi_data_2,days_kiwi,"Median transcript concentration per functional category for kiwi stages"))
annot_kiwi<-melt(t(annot_kiwi))
totalAnno_trans(annot_kiwi,kiwi_data_2,gene_mrna,"")
