library(readxl)

data<-file.choose()

if (grepl("\.xls[x]*",data,perl = T)){
  data_trans<-read_xlsx(path = data,sheet = "Transcrits")
  data_prot<-read_xlsx(path = data,sheet = "Proteines")
  colnames(data_prot)[1]<-"Protein"
  colnames(data_trans)[1]<-"Transcrit"
}

lista<-list()
total_data<-merge(data_trans,data_prot)
for (i in nrow(total_data)){
  lista[[i]]<-c("Protein_ID"=total_data[i,"Protein"],"Transcrit_ID"=total_data[i,"Transcrit"],"Transcrit_val"=total_data[i,3:29],"Protein_val"=total_data[i,30:ncol(total_data)])
}