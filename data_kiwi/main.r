setwd("../data_kiwi/")
source("../model/global.r")
# debug(solgss_Borne)
debug(confEllipse)
poids_kiwi<-read_csv("poids_kiwi.csv",col_names=c("t","y"))
# poids_kiwi<-loadData("poids_kiwi.csv","","",T))
test_data<-loadData(data = "test.xlsx",trans_sheet = "Transcrits",prot_sheet = "Proteines",F)
test_list<-test_data$parse
days_kiwi<-test_list[[1]][["DPA"]]
test_list<-Filter(function(x) x[["Transcrit_ID"]] %in% names,test_list)
# test_list<-test_list[-1]
fitWe<<-"double_sig"
coef_poids<-fitPoids(poids_kiwi[,1],poids_kiwi[,2],fitWe)
coefs_poids[["graph"]]
poids_coef<<-coef_poids$coefs
formula_poids<<-coef_poids$formula
val_mu<-mu(c(poids_kiwi$t),fitWe,poids_coef,formula_poids,dpa_analyse = NULL)
data_mu<-data.frame("DPA"=c(poids_kiwi$t),"Mu"=val_mu)
g_mu<-ggplot(data_mu,aes(x=DPA,y=Mu))+geom_line()+theme+xlab("DPA")+ylab(bquote("Growth rate "~(days^-1)))
data_rel_mu<-data.frame("DPA"=c(poids_kiwi$t),"RGR"=val_mu/fitted(coef_poids[["formula"]]))
g_rel_mu<-ggplot(data_rel_mu,aes(x=DPA,y=RGR))+geom_line()+theme+xlab("DPA")+ylab(bquote("Relative growth rate "~(days^-1)))
ggarrange(g_mu,g_rel_mu,ncol = 2)
ksmin=3*4*3*3.6*24
fitR<-"3_deg_log"
# numCores<-detectCores() - 1
# cl <- makeCluster(numCores)
# clusterEvalQ(cl, {
#   ## set up each worker.  Could also use clusterExport()
#   source("../model/global.r")
#   library(ggplot2)
#   library(grid)
#   library(egg)
#   NULL
# })
# 
# clusterExport(cl,c("poids_coef","formula_poids","ksmin","fitR","fitWe"))
# if(Sys.info()["sysname"]=="Windows"){
#   num_cor<-cl1
# }
# else if(Sys.info()["sysname"]=="Linux"){
#   num_cor<-detectCores()-1
# }
for (el in test_list){
# res_list<-pblapply(X=test_list,function(el){
  tryCatch({
    # cont<-cont+1
    print(el[["Transcrit_ID"]])
    bound_ks<-c(4.5e-3*mean(el$Transcrit_val,na.rm = T)/mean(el$Protein_val,na.rm = T),1440*mean(el$Transcrit_val,na.rm = T)/mean(el$Protein_val,na.rm = T))
    norm_data<-normaMean(el$Protein_val,el$Transcrit_val,ksmin)
    fittedmrna<<-fit_testRNA(el$DPA,norm_data$mrna,fitR)
    el$plot_mrna<-plotFitmRNA(el$DPA,norm_data$mrna,solmRNA(el$DPA,fittedmrna$coefs,fitR))
    el[["prot_mrna"]]
    par_k<-solgss_Borne(el$DPA,as.vector(norm_data$prot),as.numeric(norm_data$ks),bound_ks,"LM")
    par_k[["confEllipsePlot"]]
    if (!is.null(par_k)){
      res<-list()
      par_k[["plot_fit_prot"]]<-plotFitProt(el$DPA,as.vector(norm_data$prot),par_k$prot_fit)
      par_k[["plot_fit_prot"]]
      X<-matrice_sens(el$DPA,par_k[["solK"]][,1])
      diff<-(par_k[["error"]][["errg"]][1]*norm(as.vector(norm_data$prot),"2"))^2
      par_k[["corr_matrix"]]<-matrice_corr(X,length(norm_data$prot),diff)
      # para_min<-fminunc(par_k[["solK"]][,1],fn=minSquares,time=el$DPA,exp_data=as.vector(norm_data$prot))
      el$SOL<-par_k
      res[["TranscritID"]]<-el[["Transcrit_ID"]]
      res[["Weight formula"]]<-"Double sigmoid"
      res[["mRNA formula"]]<-fitR
      res[["Mean mRNA concentration"]]<-mean(el[["Transcrit_val"]],na.rm = T)
      res[["Mean protein concentration"]]<-mean(el[["Protein_val"]],na.rm = T)
      res[["Starting protein concentration value"]]<-unname(el[["SOL"]][["solK"]][1,1])
      res[["ks"]]<-unname(el[["SOL"]][["solK"]][2,1])
      res[["Normalized ks"]]<-unname(el[["SOL"]][["solK"]][2,1])*res[["Mean protein concentration"]]/res[["Mean mRNA concentration"]]
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
      #   write.csv(test_list[[cont]][["SOL"]][["solK"]],paste("solK/",paste(test_list[[cont]][["Transcrit_ID"]],"_Sol_ks_kd.csv"),sep = ""))
    }
    print(paste("Process finished for ",el[["Transcrit_ID"]],sep=""))
    res
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}#,cl=num_cor)

# stopCluster(cl)
valid_res<-Filter(function(x) {length(x) > 0}, res_list)
del_results<-Filter(function(x) {length(x) == 0}, res_list)
final_table<-rbindlist(valid_res)

f <- list.files("solK", include.dirs = F, full.names = T, recursive = T)
file.remove(f)
unlink("./solK", recursive = TRUE)

# save(test_list,valid_res,del_results,file=path.expand("./resultsv1.RData"))
