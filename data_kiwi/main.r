# setwd("D:/Stage M2/Stage_M2/data_kiwi/")
setwd("/media/juanma/JUANMA/Stage M2/Stage_M2/data_kiwi/")
source("../model/global2.r")

poids_kiwi<-read_csv("poids_kiwi.csv",col_names=c("t","y"))
# poids_kiwi<-loadData("poids_kiwi.csv","","",T)
days_kiwi<-rep(c(0,13,26,39,55,76,118,179,222), each = 3)
test_data<-loadData(data = "test.xlsx",trans_sheet = "Transcrits",prot_sheet = "Proteines",F)
test_list<-test_data$parse
# test_list<-sample(test_list,5)
fitWe<<-"double_sig"
coef_poids<-fitPoids(poids_kiwi[,1],poids_kiwi[,2],fitWe)
poids_coef<<-coef_poids$coefs
formula_poids<<-coef_poids$formula
val_mu<-mu(c(poids_kiwi$t),fitWe,poids_coef,formula_poids,dpa_analyse = NULL)
plot(poids_kiwi$t,val_mu,"l")
ksmin=3*4*3*3.6*24
score=0
fitR<-"3_deg_log"
cont<-0
dir.create("solK")
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
# clusterExport(cl,c("poids_coef","formula_poids","ksmin","fitR"))
if (Sys.info()["sysname"]=="Windows"){
  for (el in test_list){
  # res_list<-pblapply(X=test_list,function(el){
    tryCatch({
      # cont<-cont+1
      print(el)
      el[["Protein_val"]]<-na.omit(el[["Protein_val"]])
      print(el[["Transcrit_ID"]])
      norm_data<-normaMean(el$Protein_val,el$Transcrit_val,ksmin)
      fittedmrna<<-fit_testRNA(el$DPA,norm_data$mrna,fitR)
      el$plot_mrna<-plotFitmRNA(el$DPA,norm_data$mrna,solmRNA(el$DPA,fittedmrna$coefs,fitR))
      par_k<-solgss_Borne(el$DPA,as.vector(norm_data$prot),as.numeric(norm_data$ks),score)
      if (!is.null(par_k)){
        res<-list()
        par_k[["plot_fit_prot"]]<-plotFitProt(el$DPA,as.vector(norm_data$prot),par_k$prot_fit)
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
  }#,cl=cl)
} else if(Sys.info()["sysname"]=="Linux"){
  for (el in test_list){
  # res_list<-pblapply(X=test_list,function(el){
    # tryCatch({
    # cont<-cont+1
    print(el)
    el[["Protein_val"]]<-na.omit(el[["Protein_val"]])
    print(el[["Transcrit_ID"]])
    norm_data<-normaMean(el$Protein_val,el$Transcrit_val,ksmin)
    fittedmrna<<-fit_testRNA(el$DPA,norm_data$mrna,fitR)
    el$plot_mrna<-plotFitmRNA(el$DPA,norm_data$mrna,solmRNA(el$DPA,fittedmrna$coefs,fitR))
    par_k<-solgss_Borne(el$DPA,as.vector(norm_data$prot),as.numeric(norm_data$ks),score)
      if (!is.null(par_k)){
        el$SOL<-par_k
        res<-list()
        par_k[["plot_fit_prot"]]<-plotFitProt(el$DPA,as.vector(norm_data$prot),par_k$prot_fit)
        X<-matrice_sens(el$DPA,par_k[["solK"]][,1])
        diff<-(par_k[["error"]][["errg"]][1]*norm(as.vector(norm_data$prot),"2"))^2
        par_k[["corr_matrix"]]<-matrice_corr(X,length(norm_data$prot),diff)
        # para_min<-fminunc(par_k[["solK"]][,1],fn=minSquares,time=el$DPA,exp_data=as.vector(norm_data$prot))
        el$confEllipse<-confidenceEllipse(el[["SOL"]][["modelList"]][["model1"]],which.coef = c("ks","kd"),fill = T,segments = 52)
        if (any(el$confEllipse<0)){
          el$confEllipsePlot<-ggplot(as.data.frame(el[["confEllipse"]]),aes(x,y))+geom_path()+theme+xlim(c(0,max(as.data.frame(el[["confEllipse"]])$x)))+ylim(0,max(as.data.frame(el[["confEllipse"]])$y))+ylab("kd")+xlab("ks")
        }
        else {
          el$confEllipsePlot<-ggplot(as.data.frame(el[["confEllipse"]]),aes(x,y))+geom_path()+theme+ylab("kd")+xlab("ks")
        }
        res[["TranscritID"]]<-el[["Transcrit_ID"]]
        res[["Weight formula"]]<-"Double sigmoid"
        res[["Weight error"]]<-coef_poids[["error"]]
        res[["mRNA formula"]]<-fitR
        res[["mRNA error"]]<-fittedmrna[["error"]]
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
      # },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }#,cl=numCores)
}
# stopCluster(cl)
valid_res<-Filter(function(x) {length(x) > 0}, res_list)
del_results<-Filter(function(x) {length(x) == 0}, res_list)
final_table<-rbindlist(valid_res)

f <- list.files("solK", include.dirs = F, full.names = T, recursive = T)
file.remove(f)
unlink("./solK", recursive = TRUE)

# save(test_list,valid_res,del_results,file=path.expand("./resultsv1.RData"))
