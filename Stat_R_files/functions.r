
getData<-function(filename,sep){
  if (sep==","){
    data_fruit<-read.csv(filename,check.names = FALSE)
    data_fruit2<-data_fruit[,-1]
    rownames(data_fruit2)<-data_fruit[,1]
    delete<-rep(TRUE,nrow(data_fruit2))
    suma<-apply(data_fruit2,1,function(x) sum(x!=0))
    for (i in seq(1,length(rownames(data_fruit2)))){
      cont=1
      cont2=0
      while (cont<length(data_fruit2[i,])){
        if (all(data_fruit2[i,seq(cont,cont+2)]!=0)){
          delete[i]<-FALSE
          break
        }
        # if (sum(data_fruit2[i,cont:cont+2]!=0)==2){
        #   cont2<-cont2+1
        # }
        # if ((cont2==2) & (suma[i]!=1)){
        #   delete[i]<-FALSE
        #   break
        # }
        cont=cont+3
        
      }
    }
    index<-which(delete)
    data_fruit2<-data_fruit2[-index,]
    rownames(data_fruit2)<-data_fruit[,1][-index]
  }
  else {
    data_fruit<-read_excel(filename)
    rownames(data_fruit)<-data_fruit$...1
    data_fruit2<-data_fruit[,-1]
    rownames(data_fruit2)<-rownames(data_fruit)
    delete<-rep(TRUE,nrow(data_fruit2))
    suma<-apply(data_fruit2,1,function(x) sum(x!=0))
    for (i in seq(1,length(rownames(data_fruit2)))){
      cont=0
      pos<-c(1,3,4,6,7,9,10,12,13,15,16,18,19,20,21,23,24,26)
      for (j in seq(1,length(pos),2)){
        if (all(data_fruit2[i,seq(pos[j],pos[j+1])]!=0)){
          delete[i]<-FALSE
          break
        }
        # if (sum(data_fruit2[i,pos[j]:pos[j+1]]!=0)==2){
        #   cont<-cont+1
        # }
        # if ((cont==2) & (suma[i]!=1)){
        #   delete[i]<-FALSE
        #   break
        # }
      }
    }
    index<-which(delete)
    data_fruit2<-data_fruit2[-index,]
    rownames(data_fruit2)<-rownames(data_fruit)[-index]
  }
  # data_fruit2<-data_fruit2[apply(data_fruit2, 1, function (x) !all(x==0)),]
  data2 <- data.frame(t(data_fruit2))
  colnames(data2) <- rownames(data_fruit2)
  return(list("format1"=data_fruit2,"format2"=data2))
}
plotDensity<-function(dataframe, title,lab_x,lab_y,moyenne){
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
  if (moyenne){
    melted<-melt(colMeans(dataframe,na.rm = T))
  }
  else{
    melted<-melt(dataframe)
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

filtrage<-function(dataframe){
  dataframe[dataframe==0]<-NaN
  dataframe<-as.array(dataframe)
  index<-array()
  for (i in seq(1,length(dataframe),3)){
    if (all(is.na(dataframe[i:i+2]))){
      return(TRUE)
    }
  }
  return(index)
}


plotDensity_v2<-function(dataframe,small_data,title,lab_x,lab_y){
  melted<-melt(rowMeans(dataframe,na.rm = T))
  melted2<-melt(rowMeans(small_data,na.rm = T))
  med<-median(melted$value[melted$value>0],na.rm = T)
  med2<-median(melted2$value[melted2$value>0],na.rm = T)
  g<-ggplot()+aes(value)+geom_histogram(data=melted,aes(value),bins = 30,alpha=0.3)+scale_x_log10()+theme+
    xlab("Concentration (fmol per gFW)")+
    ylab("Number of occurences")+ggtitle(title)+geom_histogram(data = melted2,aes(value),bins=20,alpha=0.7)+geom_vline(xintercept=med,linetype="dotdash",color="red")+
    geom_vline(xintercept=med2,linetype="dotdash",color="darkblue")+annotate(geom = "text",x=med-0.1,y=(med-0.1)*2e5,label=format(round(med,2),nsmall = 2),color="red")+
    annotate(geom = "text",x=med2+0.4,y=(med2-0.1)*2e4,label=format(round(med2,2),nsmall = 2),color="darkblue")
    
  return(g)
}
normalize<-function(dataframe){
  dataframe_center<-apply(dataframe,2,function(x) x-mean(x,na.rm = T))
  dataframe<-apply(dataframe_center,2,function(x) x/sqrt(sd(x,na.rm = T)))
  # dataframe[dataframe<1e-3]<-0
  return(dataframe)
}
stageDensity<-function(dataframe,dpa,lab_x,lab_y,moyenne){
  p<-list()
  cont=1
  for (i in seq(1,length(dpa),3)){
    if (i==length(dpa)-2){
      p[[cont]]<-(plotDensity(dataframe[seq(i,i+2),],dpa[i],"","",moyenne))+theme(axis.text.x=element_text(colour="black"),axis.ticks.x=element_line(colour="black"),plot.margin= unit(c(0, 1, 0.5, 0.5), "lines"))+scale_x_continuous(trans = "log10",limits = c(1e0,1e6))
    }
    else{
      p[[cont]]<-(plotDensity(dataframe[seq(i,i+2),],dpa[i],"","",moyenne))+theme(axis.title.x = element_blank(), 
                                 axis.text.x = element_blank(), 
                                 axis.ticks.x= element_blank(), 
                                 plot.margin= unit(c(0, 1, 0.5, 0.5), "lines"))+scale_x_continuous(trans = "log10",limits = c(1e0,1e6))
    }
    cont<-cont+1
  }
  # list_plot<-do.call(rbind,p)
  # grid.draw(ggplotGrob(list_plot))
  plot_grid(plotlist = p,ncol = 1,align = "hv",hjust=1,vjust=1)
  ggarrange(plots = p,ncol=1,bottom = textGrob(lab_x,hjust=0,vjust=0),left =textGrob(lab_y,rot=90,vjust=1))
}
getHeatMap<-function(dataframe,title){
  data<-melt(cor(dataframe,use = "pairwise.complete.obs"))
  qplot(x=as.factor(Var1),y=as.factor(Var2),data = data,fill=value,geom="tile",xlab = "Sample",ylab = "Sample", main = title)+theme_classic()
  
}

clusteredHeatmap<-function(dataframe,title){
  dataframe[is.na(dataframe)]<-0
  data.dist<-vegdist(as.matrix(dataframe),method = "euclidean",na.rm = T)
  data.dist2<-vegdist(t(as.matrix(dataframe)),method = "euclidean",na.rm = T)
  scaleRYG <- colorRampPalette(c("yellow","white","blue"),
                               space = "rgb")(30)
  row.cluster<-hclust(data.dist,"ave")
  col.cluster<-hclust(data.dist2,"ave")
  log_val<-log(as.matrix(dataframe))
  log_val[which(!is.finite(log_val))]<-0
  heatmap.2(log_val,Rowv=as.dendrogram(row.cluster),
            Colv = as.dendrogram(col.cluster),
            labCol = colnames(log_val),
            col=scaleRYG,
            density.info="none",
            main = title,
            trace = "none"
  )
  # pos2<-locator()
  # pos2<-structure(list(x=c(0.27149971320082, 0.858971646016485),
  #                      y=c(0.861365598392473, 0.857450478257082)),
  #                 .Names=c("x","y"))
  # text(x=seq(pos2$x[1], pos2$x[2], len=ncol(dataframe)), y=rep(pos2$y[1],5),
  #      srt=2, xpd=TRUE, adj = 45,
  #      labels=colnames(dataframe))
}
getPCA<-function(dataframe,days_pca,title){
  # print(is.nan(dataframe))
  dataframe_center<-apply(dataframe,2,function(x) x-mean(x,na.rm = T))
  dataframe<-apply(dataframe_center,2,function(x) x/sqrt(sd(x,na.rm = T)))
  dataframe<-as.data.frame(dataframe)
  dataframe[is.na(dataframe)]<-0
  dataframe$PCA<-days_pca
  # dataframe$PCA<-as.factor(dataframe$PCA)
  data_pca<-prcomp(dataframe[,-ncol(dataframe)])
  df_out<-as.data.frame(data_pca$x)
  
  # plot(data_pca$x[,1],data_pca$x[,2])
  df_out$PCA<-dataframe$PCA
  percentage <- round(data_pca$sdev / sum(data_pca$sdev) * 100, 2)
  percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  colours<-rainbow(nrow(df_out))
  p<-ggplot(df_out,aes(x=PC1,y=PC2, color=PCA,label=PCA))
  p<-p+geom_point(size=1.2)+xlab(percentage[1])+ylab(percentage[2])+theme+scale_color_manual(values=colours,breaks=days_pca)+ggtitle(title)+geom_text_repel(size=4,force = 20)+theme(legend.position="none")
  p
  # data_pca<-PCA(dataframe[,-ncol(dataframe)],graph=FALSE)
  # plot.PCA(data_pca,choix="ind",habillage=3)
}

totalAnno_trans<-function(annot,datafruit,data_gene,title){
  annot$Var2<-NULL
  annot$value<-as.factor(annot$value)
  annot<-annot[-which(annot$value=="" | annot$Var1=="not assigned"),]
  colnames(annot)<-c("Definition","Transcrit")
  colnames(datafruit)<-paste(colnames(datafruit),"DPA",sep = " ")
  datafruit$Transcrit<-row.names(datafruit)
  a<-melt(datafruit)
  colnames(a)<-c("Transcrit","Stage","value")
  merged_data<-merge(annot,a)
  count_gene<-annot[annot$Transcrit %in% data_gene$Transcrit,]
  data2<-unique.data.frame(merged_data[,1:2])
  data2<-as.data.frame(table(data2$Definition))
  colnames(data2)<-c("Definition","Number of occurences in valid transcrits")
  data1<-as.data.frame(table(count_gene$Definition))
  colnames(data1)<-c("Definition","Number of occurences in all genes")
  A<-merge(data1,data2)
  A<-A[-which(A$Definition=="not assigned"),]
  A$Percentage<-(A$`Number of occurences in valid transcrits`/A$`Number of occurences in all genes`)*100
  # table_gene<-as.data.frame(table(count_gene$Definition))
  # table_gene<-table_gene[-which(table_gene$Var1=="not assigned"),]
  write.csv(A,"Enrichissement_transcrits_kiwi.csv",row.names = F)
  p1<-ggplot(A,aes(x=Definition))+geom_col(aes(y=`Number of occurences in all genes`))+theme(axis.text.x = element_text(angle = 0, hjust = 1,color = "black"),panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))+
    ylab("Number of genes")+xlab("")+ggtitle(title)+coord_flip()
  p2<-ggplot(A,aes(x=Definition))+geom_col(aes(y=Percentage))+theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.y=element_blank(),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))+
    ylab("Enrichment (%)")+coord_flip()+xlab("")
  grid.draw(cbind(ggplotGrob(p1),ggplotGrob(p2),size="last"))
}
totalAnno_prot<-function(annot_brute,all_corresp,datafruit,data_gene,title){
  annot_prot<-melt(t(annot_brute))
  annot_prot$Var2<-NULL
  colnames(annot_prot)<-c("Definition","Transcrit")
  annot<-annot_prot[annot_prot$Transcrit %in% all_corresp$Transcrit,]
  annot$Transcrit<-as.factor(annot$Transcrit)
  annot<-annot[-which(annot$Transcrit=="" | annot$Definition=="not assigned"),]
  # colnames(datafruit)<-paste(colnames(datafruit),"DPA",sep = " ")
  # datafruit$Transcrit<-row.names(datafruit)
  # a<-melt(datafruit)
  colnames(datafruit)<-c("Transcrit","Stage","value")
  merged_data<-merge(annot,datafruit)
  count_gene<-annot_prot[annot_prot$Transcrit %in% data_gene$Transcrit,]
  data2<-unique.data.frame(merged_data[,1:2])
  data2<-as.data.frame(table(data2[,c("Definition")]))
  colnames(data2)<-c("Definition","Number of occurences in valid proteins")
  data1<-as.data.frame(table(count_gene$Definition))
  colnames(data1)<-c("Definition","Number of occurences in all genes")
  A<-merge(data1,data2)
  # table_gene<-as.data.frame(table(count_gene$Definition))
  A$Percentage<-(A$`Number of occurences in valid proteins`/A$`Number of occurences in all genes`)*100
  A<-A[-which(A$Definition=="not assigned"),]
  write.csv(A,"Enrichissement_proteines_kiwi.csv",row.names = F)
  p1<-ggplot(A,aes(x=Definition))+geom_col(aes(y=`Number of occurences in all genes`))+theme(axis.text.x = element_text(angle = 0, hjust = 1,color = "black"),panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))+
    ylab("Number of genes")+xlab("")+ggtitle(title)+coord_flip()
  p2<-ggplot(A,aes(x=Definition))+geom_col(aes(y=Percentage))+theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.y=element_blank(),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))+
    ylab("Enrichment (%)")+coord_flip()+xlab("")
  grid.draw(cbind(ggplotGrob(p1),ggplotGrob(p2),size="last"))
}
medianAnno_prot<-function(anno_data,datafruit,days_kiwi,title){
  melt_conc<-melt(t(anno_data))
  melt_conc$Var2<-NULL
  melt_conc<-melt_conc[-which(melt_conc$value==""),]
  colnames(melt_conc)<-c("Definition","Transcrit")
  # datafruit$index<-paste(datafruit$index,"DPA",sep=" ")
  # datafruit$Transcrit<-rownames(datafruit)
  # a<-melt(datafruit)
  colnames(datafruit)<-c("Transcrit","Stage","value")
  merged_data<-merge(melt_conc,datafruit)
  merged<-aggregate(merged_data$value,by=list(merged_data$Stage,merged_data$Definition),function(x) median(x,na.rm = T))
  colnames(merged)<-c("Stage","Definition","Median")
  merged<-merged[-which(merged$Definition=="not assigned"),]
  ggplot(merged, aes(fill=factor(paste(Stage,"DPA",sep=" "),levels = days_kiwi), y=Median, x=Definition)) +
    geom_bar(position="dodge", stat="identity")+theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))+ylab("Median protein concentration (fmolgFW^1)")+
    xlab("")+ggtitle(title)+guides(fill=guide_legend(title="Stage"))
}

medianAnno_trans<-function(anno_data,datafruit,days_pca,title){
  melt_conc<-melt(t(anno_data))
  melt_conc$Var2<-NULL
  melt_conc<-melt_conc[-which(melt_conc$value==""),]
  colnames(melt_conc)<-c("Definition","Gene")
  datafruit$Gene<-row.names(datafruit)
  colnames(datafruit)<-days_pca
  a<-melt(datafruit)
  colnames(a)<-c("Gene","Stage","value")
  merged_data<-merge(melt_conc,a)
  merged<-aggregate(merged_data$value,by=list(merged_data$Stage,merged_data$Definition),median)
  colnames(merged)<-c("Stage","Definition","Median")
  p<-ggplot(merged, aes(fill=Stage, y=Median, x=Definition)) +
    geom_bar(position="dodge", stat="identity")+theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))+ylab("Median mRNA concentration (fmol·gFW^-1)")+
    xlab("")+ggtitle(title)+guides(fill=guide_legend(title="Stage"))
  print(p)
}