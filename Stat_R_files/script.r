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

# donnees_conc<- getData("Conc_Concombre.csv",",")
# 
# days_conc<-rep(c("0 DPA", "2 DPA", "5 DPA", "8 DPA", "12 DPA", "14 DPA", "18 DPA", "25 DPA", "29 DPA"), each = 3)
donnees_kiwi<-getData("Conc_Red5.csv",",")
days_kiwi<-rep(c("0 DPA", "13 DPA", "26 DPA", "39 DPA", "55 DPA", "76 DPA", "118 DPA", "179 DPA", "222 DPA"), each = 3)

# conc_data<-donnees_conc$format2
# conc_data_2<-donnees_conc$format1

kiwi_data<-donnees_kiwi$format2
kiwi_data_2<-donnees_kiwi$format1

list_tom<-getData("GSE128739_devTom_fmol.gFW.xlsx","\t")
# tom_data<-list_tom$format2
# tom_data2<-list_tom$format1
# days_tom<-c(rep(c("8 DPA", "15 DPA", "21 DPA", "28 DPA", "34 DPA", "41 DPA"),each=3),"49 DPA","49 DPA",rep(c("50 DPA","53 DPA"),each=3))
rm(donnees_conc,donnees_kiwi,list_tom)
# my_plots <- lapply(rownames(kiwi_data)[1116:1125], function(var_x){
#   p <-
#     ggplot(kiwi_data)+
#     aes_string(var_x)+
#     geom_line(kiwi_data[c(var_x),])
#
# })
# plot_grid(plotlist=my_plots)

#ggplot(conc_data.melt, aes(x = variable, y = Gene, fill =value)) + geom_tile()
# conc_data.melt<-melt(conc_data,id.vars = "Gene")
#
# ggplot(conc_data.melt, aes(x=value)) + geom_bar() + xlim(0,10)

#conc_data<-melt(conc_data)
#ggplot(conc_data, aes(variable, value, group=factor(Gene))) + geom_line(aes(color=factor(Gene)))


# pc<-prcomp(kiwi_data)
# names(pc)
# plot(pc$x[, 1], pc$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

# Densities ---------------------------------------------------------------
# tom_data[tom_data<1e-10]<-0
# stageDensity(conc_data,days_conc,"Distribution of cucumber concentration")
stageDensity(kiwi_data,days_kiwi,"Concentration (fmol per gFW)","Number of variables",moyenne = T)
# pos<-c(1,3,4,6,7,9,10,12,13,15,16,18,19,20,21,23,24,26)
# cont<-1
# for (j in seq(1,length(pos),2)){
#   title<-paste("Distribution of tomato concentration",unique(days_tom)[cont],sep = ": ")
#   print(plotDensity(tom_data[pos[j]:pos[j+1],],title))
#   cont<-cont+1
# }
#

# for (i in seq(1,length(days_conc),3)){
#   var_name<-gsub(" ","_",days_conc[i])
#   assign(sprintf("gplot_%s",var_name),plotDensity(conc_data[i:i+2,],title))
# }
# plots<-ls(pattern = "gplot_*")
# print(plots)
# multiplot(gplot_0_DPA,gplot_2_DPA,gplot_5_DPA,gplot_8_DPA,gplot_12_DPA,gplot_14_DPA,gplot_18_DPA,gplot_25_DPA,gplot_29_DPA)

# g1<-plotDensity(conc_data,"Distribution of cucumber concentration")
g2<-plotDensity(kiwi_data,"","Concentration (fmol per gFW)","Number of occurences",T)
# g3<-plotDensity(tom_data,"Distribution of tomato concentration")

# print(g1)
print(g2)
# print(g3)
# 
# melted_conc<-melt(colMeans(conc_data))
# melted_kiwi<-melt(colMeans(kiwi_data))
# melted_tom<-melt(colMeans(tom_data))
# 
# g4<-ggplot()+
#   geom_histogram(melted_conc,mapping=aes(x=value,fill="BAR YELLOW"),bins = 50,alpha=0.5,fill="yellow")+scale_x_log10()+xlab("Concentration (fmol per gFW)")+
#   ylab("Number of variables") +theme+
#   geom_vline(xintercept=median(melted_conc$value[melted_conc$value>0]),linetype="dotdash",color="yellow")+
#   geom_histogram(melted_kiwi,mapping=aes(x=value,fill="BAR GREEN"),bins = 50,alpha=0.2,fill="green")+
#   geom_vline(xintercept=median(melted_kiwi$value[melted_kiwi$value>0]),linetype="dotdash",color="green")+
#   geom_histogram(melted_tom,mapping=aes(x=value,fill="BAR RED"),bins = 50,alpha=0.3,fill="red")+scale_fill_manual("Species",values=c("red","green","blue"),labels=c("Cucumber","Kiwi","Tomato"))+
#   geom_vline(xintercept=median(melted_tom$value[melted_tom$value>0]),linetype="dotdash",color="red")+ggtitle("Distribution of concentration for kiwi, cucumber and tomato")
# print(g4)
# Heatmap -----------------------------------------------------------------
# heatmaply(as.matrix(conc_data_2))

getHeatMap(conc_data_2,"Correlation heatmap for cucumber")
getHeatMap(kiwi_data_2,"Correlation heatmap for kiwi Red5")
# getHeatMap(tom_data2,"Correlation heatmap for tomato")
#
clusteredHeatmap(conc_data_2,"Heatmap of conc for cucumber")
clusteredHeatmap(kiwi_data_2,"Heatmap of conc for kiwi")
    # clusteredHeatmap(tom_data2,"Heatmap of conc for tomato")
# #---------------------------------------------------------------------

getPCA(conc_data,days_conc,"ACP of cucumber concentration")
getPCA(kiwi_data,days_kiwi,"ACP of kiwi concentration")
# getPCA(tom_data,days_tom,"ACP of tomato concentration")
# #
# stage_conc<-aggregate(conc_data,list(rep(1:(nrow(conc_data)%/%3+1),each=3,len=nrow(conc_data))),mean)[-1]
# rownames(stage_conc)<-unique(days_conc)
stage_kiwi<-aggregate(kiwi_data,list(rep(1:(nrow(kiwi_data)%/%3+1),each=3,len=nrow(kiwi_data))),function(x) mean(x,na.rm = T))[-1]
rownames(stage_kiwi)<-unique(days_kiwi)

# pos<-c(1,3,4,6,7,9,10,12,13,15,16,18,19,20,21,23,24,26)
# stage_tom<-data.frame()
# for (j in seq(1,length(pos),2)){
#   data<-apply(tom_data[pos[j]:pos[j+1],], 2, mean)
#   stage_tom<-rbind(stage_tom,data)
# }
# rownames(stage_tom)<-unique(days_tom)
# #
getPCA(stage_conc,unique(days_conc),"ACP of cucumber concentration: average of samples for each stade")
getPCA(stage_kiwi,unique(days_kiwi),"PCA of mRNA concentration")
  # getPCA(stage_tom,unique(days_tom),"ACP of tomato concentration: average of samples for each stade")

varr_conc<-apply(stage_conc,1,function(x) cv(x)*100)
var_kiwi<-apply(stage_kiwi,1,function(x) cv(x)*100)
var_tom<-apply(stage_tom,1,function(x) cv(x)*100)
data<-data.frame(varr_conc,var_kiwi,var_tom)
colnames(data)<-c("Concombre","Kiwi","Tomate")
data<-melt(data)
data$stade<-as.factor(rep(seq(1:9),3))
colnames(data)<-c("Espèce","value","stade")
ggplot(data, aes(fill=Esp?ce, y=value, x=stade)) +
  geom_bar(position="dodge", stat="identity")+theme+scale_fill_manual(values=c("yellow","green","red"))+ylab("CV (%)")



# Annot par stade ---------------------------------------------------------

annot_kiwi<-read.csv("Prot_red5.csv",header = T,sep = "\t",check.names = F)
annot_kiwi<-annot_kiwi[,-1]
annot_conc<-read.csv("Annot_genes_conc.txt",header = T,sep = "\t",check.names = F)
annot_conc<-annot_conc[,-1]

gene_mrna<-read.csv("gene_mrna_Red5.csv",check.names = F)

# pdf("Medianne_annotations_Red5.pdf",width=8.5,height = 11)
print(medianAnno_trans(annot_kiwi,kiwi_data_2,days_kiwi,"Median transcript concentration per functional category for kiwi stages"))
annot_kiwi<-melt(t(annot_kiwi))
totalAnno_trans(annot_kiwi,kiwi_data_2,gene_mrna,"")
# dev.off()