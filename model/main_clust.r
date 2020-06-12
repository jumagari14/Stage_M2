
# setwd("D:/Stage M2/Stage_M2/model/")
# setwd("/media/juanma/JUANMA/Stage M2/Stage_M2/model/")
library(getopt
)
spec <- matrix(c(
  "workDir", "o", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)
setwd(opt$workDir)
source("functions.r")

# test_data<-lista[[30]]

poids<-read_csv("poids_test.csv",col_names=c("DPA","Poids"))
poids_kiwi<-read_xlsx("Kiwi_FW.xlsx",sheet = "Kiwifruit")
per_dpa<-days_kiwi<-rep(c(0,13,26,39,55,76,118,179,222), each = 3)
test_data<-loadData(data = "test.xlsx",trans_sheet = "Transcrits",prot_sheet = "Proteines",F)
test_list<-test_data$parse
test_list<-sample(test_list,3)
coef_poids<-fitPoids(poids_kiwi$DPA,poids_kiwi$Weight_g,"double_sig")
poids_coef<<-coef_poids$coefs
formula_poids<<-coef_poids$formula
val_mu<-mu(c(poids_kiwi$DPA),"double_sig",poids_coef,formula_poids,dpa_analyse = NULL)
plot(poids_kiwi$DPA,val_mu,"l")
ksmin=3*4*3*3.6*24
score=0
cont<<-0
dir.create("solK")
numCores <- detectCores()-1
# cl <- makeCluster(detectCores()-1, type='PSOCK')
# registerDoParallel(numCores)
res_list<-mclapply(test_list,function(el){
  el[["Protein_val"]]<-na.omit(el[["Protein_val"]])
  print(cont)
  norm_data<-normaMean(el$Protein_val,el$Transcrit_val,ksmin)
  fitR<<-"3_deg_log"
  fittedmrna<<-fit_testRNA(el$DPA,norm_data$mrna,fitR)
  el$plot_mrna<-plotFitmRNA(el$DPA,norm_data$mrna,solmRNA(el$DPA,fittedmrna,"3_deg"))
  par_k<-solgss_Borne(el$DPA,as.vector(norm_data$prot),as.numeric(norm_data$ks),score)
  if (!is.null(par_k)){
    par_k[["plot_fit_prot"]]<-plotFitProt(el$DPA,as.vector(norm_data$prot),par_k$prot_fit)
    X<-matrice_sens(el$DPA,par_k[["solK"]][,1])
    diff<-(par_k[["error"]][["errg"]][1]*norm(as.vector(norm_data$prot),"2"))^2
    par_k[["corr_matrix"]]<-matrice_corr(X,length(norm_data$prot),diff)
    para_min<-fminunc(par_k[["solK"]][,1],fn=minSquares,time=el$DPA,exp_data=as.vector(norm_data$prot))
    el$SOL<-par_k
    if (any(para_min$par<0)){
      init_prot<-init_conc(el$DPA,as.vector(norm_data$prot))
      para_min<-fmincon(par_k[["solK"]][,1],fn=minSquares,time=el$DPA,exp_data=as.vector(norm_data$prot),lb=c(init_prot$min,0,0),ub=c(init_prot$max,Inf,Inf))
    }
    #   write.csv(test_list[[cont]][["SOL"]][["solK"]],paste("solK/",paste(test_list[[cont]][["Transcrit_ID"]],"_Sol_ks_kd.csv"),sep = ""))
  }
  el
  # },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
},mc.cores = numCores)

valid_res<-Filter(function(x) {length(x) > 6}, res_list)
del_results<-Filter(function(x) {length(x) == 6}, res_list)
files2zip <- dir('solK', full.names = TRUE)
zip(zipfile = 'testZip', files = files2zip)

f <- list.files("solK", include.dirs = F, full.names = T, recursive = T)
file.remove(f)
unlink("./solK", recursive = TRUE)

# save(test_list,valid_res,del_results,file=path.expand("./resultsv1.RData"))
