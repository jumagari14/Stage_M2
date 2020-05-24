library(readxl)
setwd("/media/juanma/JUANMA/Stage M2/Stage_M2/model/")

# data<-file.choose()
loadData<-function(data,trans_sheet,prot_sheet){
  if (grepl("xls[x]*",data,perl = T)){
    data_trans<-read_xlsx(path = data,sheet = trans_sheet)
    data_prot<-read_xlsx(path = data,sheet =prot_sheet)
    colnames(data_prot)[1]<-"Protein"
    colnames(data_trans)[1]<-"Transcrit"
  }
  
  t<-round(as.numeric(colnames(data_prot)[-which(is.na(as.numeric(as.character(colnames(data_prot)))))]))
  
  total_data<-merge(data_trans,data_prot)
  lista<-vector("list",nrow(data_trans))
  for (i in seq(1,nrow(total_data))){
    lista[[i]]<-list("Protein_ID"=total_data[i,"Protein"],"Transcrit_ID"=total_data[i,"Transcrit"],"Transcrit_val"=as.matrix(total_data[i,3:29]),"Protein_val"=as.matrix(total_data[i,30:ncol(total_data)]),"DPA"=t)
    
  }
  return(list("prot"=data_prot,"mrna"=data_trans,"parse"=lista))
}