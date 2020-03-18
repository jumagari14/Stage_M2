
# Load packages -----------------------------------------------------------
#setwd("E:/Stage M2/Analyse_stats")
setwd("/media/juanma/JUANMA/Stage M2/Analyse_stats/")
require("readxl")
require("gplots")
require("ggfortify")
require("GGally")
require("FactoMineR")
require("factoextra")
require("reshape2")
require("cowplot")
require("dplyr") 
require("vegan")
require("stringr") 
require("tidyr")
require("plotly")
require("heatmaply")

# Load Data ---------------------------------------------------------------


getData<-function(filename,sep){
  if (sep==","){
    data_fruit<-read.csv(filename)
    data_fruit2<-data_fruit[,-1]
    rownames(data_fruit2)<-data_fruit[,1]
  }
  else {
    data_fruit<-read_excel(filename)
    data_fruit2<-data_fruit[,-1]
    rownames(data_fruit)<-data_fruit$..1
  }
  data_fruit2<-data_fruit2[apply(data_fruit2[,-1], 1, function(x) !all(x==0)),]
  data2 <- data.frame(t(data_fruit2))
  colnames(data2) <- rownames(data_fruit2)
  return(list("format1"=data_fruit2,"format2"=data2))
}

donnees_conc<- getData("Conc_Concombre.csv",",")

days_conc<-rep(c("0 DPA", "2 DPA", "5 DPA", "8 DPA", "12 DPA", "14 DPA", "18 DPA", "25 DPA", "29 DPA"), each = 3)
donnees_kiwi<-getData("Conc_Kiwifruit.csv",",")
days_kiwi<-rep(c("0 DPA", "13 DPA", "26 DPA", "39 DPA", "55 DPA", "76 DPA", "118 DPA", "179 DPA", "222 DPA"), each = 3)

conc_data<-donnees_conc$format2
kiwi_data<-donnees_kiwi$format2

conc_data_2<-donnees_conc$format1
kiwi_data_2<-donnees_kiwi$format1
list_tom<-getData("GSE128739_devTom_fmol.gFW.xlsx","\t")
tom_data<-list_tom$format2
tom_data2<-list_tom$format1
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

plotDensity<-function(dataframe){
  melted<-melt(dataframe)
  ggplot(melted)+aes(value)+geom_histogram(bins = 30)+scale_x_log10()+xlab("Concentration (fmol per gFW)")+ 
    ylab("Number of variables") +theme+ #+ 
    # stat_summary(aes(x=value,group=1),fun.y = median,
    #             geom = "line", width = 0.5,group=1)
    geom_vline(xintercept=median(melted$value[melted$value!=0]),linetype="dotdash",color="red")
}

plotDensity(conc_data)
plotDensity(kiwi_data)
plotDensity(tom_data)

# Heatmap -----------------------------------------------------------------

getHeatMap<-function(dataframe){
  data<-melt(cor(dataframe))
  qplot(x=Var1,y=Var2,data = data,fill=value,geom="tile")+theme_classic()
  
}

clusteredHeatmap<-function(dataframe){
  data.dist<-vegdist(as.matrix(dataframe),method = "euclidean")
  scaleRYG <- colorRampPalette(c("yellow","white","blue"), 
                               space = "rgb")(30) 
  row.cluster<-hclust(data.dist,"aver")
  log_val<-log(as.matrix(dataframe))
  log_val[which(!is.finite(log_val))]<-0
  heatmap.2(log_val,Rowv=as.dendrogram(row.cluster),
            Colv = F,
            labRow = F,
            col=scaleRYG,
            margin=c(5,5),
            density.info="none",
            trace = "none"
  )
  # pos2<-locator()
  pos2<-structure(list(x=c(0.27149971320082, 0.858971646016485),
                       y=c(0.861365598392473, 0.857450478257082)),
                  .Names=c("x","y"))
  text(x=seq(pos2$x[1], pos2$x[2], len=ncol(dataframe)), y=rep(pos2$y[1],5),
       srt=2, xpd=TRUE, adj = 45,
       labels=colnames(dataframe))
}
# heatmaply(as.matrix(conc_data_2))

getHeatMap(conc_data_2)
getHeatMap(kiwi_data_2)
getHeatMap(tom_data2)

png("heatmap_concombre.png", 1000, 1000,pointsize = 20)
clusteredHeatmap(conc_data_2)
dev.off()
png("heatmap_kiwi.png",100,100)
clusteredHeatmap(kiwi_data_2)
dev.off()
# PCA ---------------------------------------------------------------------

getPCA<-function(dataframe,days_pca){
  dataframe$PCA=days_pca
  dataframe$PCA<-as.factor(dataframe$PCA)
  # data_pca<-prcomp(dataframe)
  # df_out<-as.data.frame(data_pca$x)
  # 
  # # plot(data_pca$x[,1],data_pca$x[,2])
  # 
  # percentage <- round(data_pca$sdev / sum(data_pca$sdev) * 100, 2)
  # percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  # p<-ggplot(df_out,aes(x=PC1,y=PC2, label=row.names(dataframe) ))
  # p<-p+geom_point()+ geom_text(size=2)+xlab(percentage[1])+ylab(percentage[2])+theme
  # p
  data_pca<-PCA(dataframe[,-ncol(dataframe)],graph=FALSE)
  plot.PCA(data_pca,choix="ind",habillage=3)
}

getPCA(conc_data,days_conc)
getPCA(kiwi_data,days_kiwi)
