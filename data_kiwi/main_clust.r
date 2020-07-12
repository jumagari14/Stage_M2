library(getopt,quietly=TRUE, warn.conflicts=FALSE)
args <- commandArgs(trailingOnly=TRUE)
spec <- matrix(c(
  "workDir", "o", 1, "character",
  "mainFile","f",1, "character",
  "weightFile","w",1, "character",
  "finalFile","a",1,"character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)
print(opt$workDir)
setwd(opt$workDir)
source("../model/global2.r")

# test_data<-lista[[30]]
pdf(NULL)
poids_kiwi<-read_xlsx(opt$weightFile,sheet = "Poids")
test_data<-loadData(data = opt$mainFile,trans_sheet = "RNA",prot_sheet = "Protein",F)
test_list<-test_data$parse
days_kiwi<-test_list[[1]][["DPA"]]
fitWe<<-"double_sig"
coef_poids<-fitPoids(poids_kiwi[,1],poids_kiwi[,2],fitWe)
poids_coef<<-coef_poids$coefs
formula_poids<<-coef_poids$formula
val_mu<-mu(c(poids_kiwi$DPA),fitWe,poids_coef,formula_poids,dpa_analyse = NULL)
plot(poids_kiwi$DPA,val_mu,"l")
ksmin=3*4*3*3.6*24
score=0
cont<-0
dir.create("solK")
numCores <- detectCores()-1
print(numCores)
# cl <- makeCluster(detectCores()-1, type='PSOCK')
# registerDoParallel(numCores)
res_list<-mclapply(test_list,function(el){
  tryCatch({
  cont<-cont+1
  print(el[["Transcrit_ID"]])
  res<-list()
  el[["Protein_val"]]<-na.omit(el[["Protein_val"]])
    norm_data<-normaMean(el$Protein_val,el$Transcrit_val,ksmin)
    fitR<<-"3_deg_log"
    fittedmrna<<-fit_testRNA(el$DPA,norm_data$mrna,fitR)
    el$plot_mrna<-plotFitmRNA(el$DPA,norm_data$mrna,solmRNA(el$DPA,fittedmrna$coefs,fitR))
    par_k<-solgss_Borne(el$DPA,as.vector(norm_data$prot),as.numeric(norm_data$ks),score)
    if (!is.null(par_k)){
      par_k[["plot_fit_prot"]]<-plotFitProt(el$DPA,as.vector(norm_data$prot),par_k$prot_fit)
      X<-matrice_sens(el$DPA,par_k[["solK"]][,1])
      diff<-(par_k[["error"]][["errg"]][1]*norm(as.vector(norm_data$prot),"2"))^2
      par_k[["corr_matrix"]]<-matrice_corr(X,length(norm_data$prot),diff)
      # para_min<-fminunc(par_k[["solK"]][,1],fn=minSquares,time=el$DPA,exp_data=as.vector(norm_data$prot))
      el$SOL<-par_k
      res[["TranscritID"]]<-el[["Transcrit_ID"]]
      res[["Weight formula"]]<-"Double sigmoid"
      res[["Weight error"]]<-coef_poids[["error"]]
      res[["mRNA formula"]]<-fitR
      res[["mRNA error"]]<-fittedmrna[["error"]]
      res[["Mean mRNA concentration"]]<-mean(el[["Transcrit_val"]],na.rm = T)
      res[["Mean protein concentration"]]<-mean(el[["Protein_val"]],na.rm = T)
      res[["Starting protein concentration value"]]<-unname(el[["SOL"]][["solK"]][1,1])
      res[["ks"]]<-unname(el[["SOL"]][["solK"]][2,1])*res[["Mean protein concentration"]]/res[["Mean mRNA concentration"]]
      res[["Normalized ks"]]<-unname(el[["SOL"]][["solK"]][2,1])
      res[["kd"]]<-unname(el[["SOL"]][["solK"]][3,1])
      res[["Fitting error value"]]<-el[["SOL"]][["error"]][["errg"]][1]
      res[["Fitting error score"]]<-el[["SOL"]][["error"]][["score"]]
      res[["Fitting error message"]]<-el[["SOL"]][["error"]][["message"]]
      res[["Optimization error score"]]<-el[["SOL"]][["opt_eval"]][["score"]]
      res[["Optimization error message"]]<-el[["SOL"]][["opt_eval"]][["message"]]
      # if (any(para_min$par<0)){
      #   init_prot<-init_conc(el$DPA,as.vector(norm_data$prot))
      #   para_min<-fmincon(par_k[["solK"]][,1],fn=minSquares,time=el$DPA,exp_data=as.vector(norm_data$prot),lb=c(init_prot$min,0,0),ub=c(init_prot$max,Inf,Inf))
      # }
      # print(res)
    #   write.csv(el[["SOL"]][["solK"]],paste(el[["Transcrit_ID"]],"_Sol_ks_kd.csv",sep = ""))
    }
  res
  },error=function(e){cat("ERROR :",conditionMessage(e)," for ",el[["Transcrit_ID"]], "\n")})
},mc.cores = numCores,mc.preschedule=TRUE)
valid_res<-Filter(function(x) {length(x) > 1}, res_list)
del_results<-Filter(function(x) {length(x) == 0}, res_list)
final_table<-rbindlist(valid_res)
dev.off()
write.csv(final_table,opt$finalFile)

# save(test_list,valid_res,del_results,file=path.expand("./resultsv1.RData"))
