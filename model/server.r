library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(grid)
library(egg)

# source("input.r")


function(input, output, session) {
  js$disableTab("tabRes")
  fit_op<-reactiveValues(data=NULL)
  run_calc<-reactiveValues(data=NULL)
  en_but<-reactiveValues(enable=FALSE)
  theme<<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
  
  output$mult_files<-renderUI({
    if (!input$multiple){
      tagList(fileInput("prot_file","Choose protein file"),
              fileInput("mrna_file","Choose transcript file"))
    }
  })
  
  output$sing_file<-renderUI({
    if (input$multiple){
      tagList(textInput("protein_tab","Name of protein tab",value = "Proteines"),
              textInput("rna_tab","Name of mRNA tab",value = "Transcrits"),
              fileInput("data_file","Choose xls/xlsx file",accept=c(".xls",".xlsx")))
    }
  })
  observeEvent(input$method_we, {
    # updateTabsetPanel(session, "params", selected = input$method_we)
    updateTabsetPanel(session,"formulas",selected = input$method_we)
  })
  observe({
    if(!is.null(input$data_file)){
      inFile<-input$data_file
      list_data<-loadData(inFile$datapath,input$rna_tab,input$protein_tab,poids=F)
      mrna_data<-list_data$mrna
      prot_data<-list_data$prot
      test_list<-list_data$parse
      test_list<<-sample(test_list,3)
      clean_mrna_data<<-mrna_data[,-which(is.na(as.numeric(as.character(colnames(mrna_data)))))]
      clean_prot_data<<-prot_data[,-which(is.na(as.numeric(as.character(colnames(prot_data)))))]
    }
  }) 
  observe({
    if((!is.null(input$prot_file)) & (!is.null(input$mrna_file))){
      protFile<-input$prot_file
      mrnaFile<-input$mrna_file
      prot_data<-loadData(protFile$datapath,"","",poids=F)
      mrna_data<-loadData(mrnaFile$datapath,"","",poids=F)
      clean_mrna_data<<-mrna_data[,-which(is.na(as.numeric(as.character(colnames(mrna_data)))))]
      clean_prot_data<<-prot_data[,-which(is.na(as.numeric(as.character(colnames(prot_data)))))]
      
      total_data<-merge(mrna_data,prot_data)
      lista<-vector("list",nrow(mrna_data))
      for (i in seq(1,nrow(total_data))){
        lista[[i]]<-list("Protein_ID"=total_data[i,"Protein"],"Transcrit_ID"=total_data[i,"Transcrit"],"Transcrit_val"=as.matrix(total_data[i,3:29]),"Protein_val"=as.matrix(total_data[i,30:ncol(total_data)]),"DPA"=t)
        
      }
      # test_list<<-lista
      test_list<<-sample(lista,3)
    }
  })    
  observeEvent(input$disp_distr,{
    print("Plotting...")
    output$distr_plot<-renderPlot({print(combineGraphs(clean_mrna_data,clean_prot_data,"",moyenne = T))})
    print("Finished")
  })
  
  # parList<-reactiveValues()
  # observe({
  #   for (i in reactiveValuesToList(input)){
  #   print(i)
  #   if (grepl("par[1-9]+_sig",i,perl = T)){
  #     newlist[[input[[i]]]]<-input[[i]]
  #   }
  # }
  # # })
  
  observe({
    parList<<-callModule(paramList,"params",input$method_we)})
  observeEvent(input$fit_op,{
    print(parList())
    inFile<-input$weight_data
    days_kiwi<-rep(c(0,13,26,39,55,76,118,179,222), each = 3)
    poids_data<-loadData(inFile$datapath,"","",poids=T)
    print("Fitting...")
    tryCatch({
      coefs_poids<<-fitPoids_v2(poids_data[,1],poids_data[,2],input$method_we,parList())
    },
    warning = function(warn){
      showNotification(paste0(warn), type = 'warning')
    },
    error = function(err){
      showNotification(paste0(err), type = 'err')
    })
    print(coefs_poids$coefs)
    val_mu<-mu(c(poids_data$DPA),input$method_we,coefs_poids$coefs,coefs_poids$formula,dpa_analyse = NULL)
    data_mu<-data.frame("DPA"=c(poids_data$DPA),"Mu"=val_mu)
    g_mu<<-ggplot(data_mu,aes(x=DPA,y=Mu))+geom_line()+theme+xlab("DPA")+ylab("Growth rate (days^-1)")
    fit_op$state<-TRUE
    print("Finished!!")
    output$fitplot<-renderPlot({
      req((fit_op$state)==TRUE,exists("coefs_poids"))
      ggarrange(coefs_poids$graph,g_mu,ncol=2)
    })
  })
  observeEvent(input$run_loop,{
    if (input$fit_mrna!=""){
      ksmin=as.numeric(as.character(input$ksmin))
      score=0
      cont<-0
      poids_coef<<-coefs_poids$coefs
      formula_poids<<-coefs_poids$formula
      mess<-showNotification(paste("Running..."),duration = NULL,type = "message")
      for (el in test_list){
        tryCatch({
          run_calc$run<-TRUE
          cont<-cont+1
          print(cont)
          norm_data<-normaMean(el$Protein_val,el$Transcrit_val,ksmin)
          fittedmrna<<-fit_testRNA(el$DPA,norm_data$mrna,"3_deg")
          par_k<-solgss_Borne(el$DPA,as.vector(norm_data$prot),as.numeric(norm_data$ks),score)
          par_k[["plot_fit_prot"]]<-plotFitProt(el$DPA,as.vector(norm_data$prot),par_k$prot_fit)
          X<-matrice_sens(el$DPA,par_k[["solK"]][,1])
          diff<-(par_k[["error"]][["errg"]][1]*norm(as.vector(norm_data$prot),"2"))^2
          par_k[["corr_matrix"]]<-matrice_corr(X,length(norm_data$prot),diff)
          if (!is.null(par_k)){
            test_list[[cont]]$SOL<-par_k
            # write.csv(test_list[[cont]][["SOL"]][["solK"]],paste("solK/",paste(test_list[[cont]][["Transcrit_ID"]],"_Sol_ks_kd.csv"),sep = ""))
          }
        },error=function(e){showNotification(paste0("Protein fitting not achieved for ",el$Transcrit_ID,sep=" "),type = "error",duration = NULL)})
        
      }
      valid_res<<-Filter(function(x) {length(x) > 5}, test_list)
      print(valid_res[[1]])
      mess<-showNotification(paste("Finished!!"),duration = NULL,type = "message")
      en_but$enable<-TRUE
      
    }
  })
  output$downFile<-downloadHandler(
    filename = function(){
      paste("results_KsKd-",Sys.Date(),".zip",sep="")
    },
    content = function(file){
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;
      for (res in valid_res){
        fileName<-paste(res[["Transcrit_ID"]],"_Sol_ks_kd.csv",sep = "")
        write.csv(res[["SOL"]][["solK"]],fileName)
        files<-c(files,fileName)
      }
      zip(zipfile = file,files = files)
      if(file.exists(paste0(file, ".zip"))) {file.rename(paste0(file, ".zip"), file)}
    },contentType = "application/zip"
    
  )
  observe({
    if (en_but$enable){
      enable("downFile")
    }
  })
}