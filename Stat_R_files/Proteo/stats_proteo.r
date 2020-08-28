# Load packages and data -----------------------------------------------------------
setwd("E:/Stage M2/Analyse_stats/Proteo/")

# setwd("/media/juanma/JUANMA/Stage M2/Analyse_stats/Proteo")
load("kiwi_juan_Conc(fmolgFW)_index_correc_sansUPS_nouvelle_tableMW.RData")
library(readxl)
library(gplots)
library(reshape2)
library(vegan)
library(ggplot2)
library(grid)
library(staplr)
library(cowplot)
library(egg)
library(ggrepel)
source("../functions.r")




# Main --------------------------------------------------------------------


calcul=c("ibaq_fmolgFW","modelssRep_fmolgFW","modelAcRep_fmolgFW")

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
days_kiwi<-rep(c("0 DPA", "13 DPA", "26 DPA", "39 DPA", "55 DPA", "76 DPA", "118 DPA", "179 DPA", "222 DPA"), each = 3)
prot_conc<-conc.sans.sdRT[conc.sans.sdRT$variable==calcul[3],]
# prot_conc[is.na(prot_conc$value),"value"]<-0
dataprot<-dcast(prot_conc,protein~index,fun.aggregate = function(x) mean(x,na.rm = T), value.var = "value")
rownames(dataprot)<-dataprot$protein
dataprot<-dataprot[,-1]
dataprot<-dataprot[,sort(as.numeric(colnames(dataprot)),index=T)$ix]
colnames(dataprot)<-unique(days_kiwi)
dataprot2<-as.data.frame(t(dataprot))
unique_prot<-unique(prot_conc$protein)
prot_conc2<-prot_conc

for (i in seq(length(unique_prot))){
  prot_sample<-prot_conc[prot_conc$protein==unique_prot[i],]
  prot_sample$index<-make.unique(prot_sample$index,sep=".")
  prot_conc2[prot_conc2$protein==unique_prot[i],]<-prot_sample
}

data_csv<-dcast(prot_conc2,protein~index,value.var = "value")
rownames(data_csv)<-data_csv$protein
data_prot_kiwi2<-data_csv[,-1]
data_prot_kiwi2<-data_prot_kiwi2[,sort(as.numeric(colnames(data_csv)),index=T)$ix]
data_prot_kiwi<-as.data.frame(t(data_prot_kiwi2))

# Density ---------------------------------------------------------------
# pdf("stats_proteo_sansUPS_nouvMW.pdf",width = 8.5,height = 11)
plotDensity(data_prot_kiwi,"","Concentration (fmol per gFW)","Number of occurences",moyenne = T)
stageDensity(data_prot_kiwi,days_kiwi,"Concentration (fmol per gFW)","Number of variables",moyenne = T)
# Heatmap -----------------------------------------------------------------

getHeatMap(dataprot,"Correlation heatmap for kiwi protein concentration")
clusteredHeatmap(dataprot,"Heatmap of protein concentration")


# PCA ---------------------------------------------------------------------

getPCA(dataprot2,unique(days_kiwi),"PCA of protein concentration")
  # dev.off()

# Annot par stade  --------------------------------------------------------
# pdf("Annot_fonctionnelle_prot.pdf")
annot_kiwi<-read.csv("../Prot_red5.csv",check.names = F,sep = "\t")
annot_kiwi<-annot_kiwi[,-1]

gene_mrna<-read.csv("../gene_mrna_Red5.csv",check.names = F)

gene_mrna_prot<-read.csv("../gene_prot_mrna_red5.csv",check.names = F)
colnames(gene_mrna_prot)[3]<-"protein"

trans_id<-merge(prot_conc[,c("protein","index","value")],gene_mrna_prot[,c("Transcrit","protein")],by="protein")
trans_prot<-trans_id[,c("Transcrit","index","value")]


print(medianAnno_prot(annot_kiwi,trans_prot,unique(days_kiwi),""))
# dev.off()
# pdf("Prueba2.pdf",,width = 8.5,height = 11)
totalAnno_prot(annot_kiwi,gene_mrna_prot,trans_prot,gene_mrna,"")
# dev.off()
# staple_pdf(input_files = c("stats_proteo_sansUPS_nouvMW.pdf","Prueba2.pdf"),output_filepath = "stats_proteo_sansUPS_nouvelleMW.pdf",overwrite = T)
# file.remove("Prueba2.pdf")
#  Save prot conc data as correct format -------------------------------------------------------------------------



colnames(data_csv)<-sub("\.","_",colnames(data_csv))
write.csv(data_csv,"Prot_conc_modelAc_sansUPS_nouvMW.csv",sep = ",")
