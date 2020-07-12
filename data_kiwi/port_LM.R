setwd("../data_kiwi/")

data_kskd<-read.csv("../data_kiwi/ver23JuinComplete_devCorrected.csv",check.names = F)
data2<-data_kskd[which((data_kskd$`Optimization error score`==10) & (data_kskd$`Fitting error score`>6)),]

data3<-read.csv("ver28JuinResults_AlgPort.csv",check.names = F)
data4<-data3[which((data3$`Optimization error score`==10) & (data3$`Fitting error score`>6)),]

data5<-merge(data2,data4,by="TranscritID")
