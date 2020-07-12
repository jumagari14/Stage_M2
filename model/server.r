library(shiny)
library(shinyjs)
library(ggplot2)
library(grid)
library(egg)

# source("input.r")


function(input, output, session) {
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
      test_el<<-sample(test_list,1)[[1]]
    }
    if((!is.null(input$prot_file)) & (!is.null(input$mrna_file))){
      protFile<-input$prot_file
      mrnaFile<-input$mrna_file
      mrna_data<-loadData(mrnaFile$datapath,"","",poids=F)
      prot_data<-loadData(protFile$datapath,"","",poids=F)
      colnames(prot_data)[1]<-"Protein"
      colnames(mrna_data)[1]<-"Transcrit"
      ind_num_mrna<-which(!is.na(as.numeric(as.character(colnames(mrna_data)))))
      ind_prot_num<-which(!is.na(as.numeric(as.character(colnames(prot_data)))))
      mrna_data[,ind_num_mrna]<-as.data.frame(sapply(mrna_data[,ind_num_mrna],as.numeric))
      prot_data[,ind_prot_num]<-as.data.frame(sapply(prot_data[,ind_prot_num],as.numeric))
      clean_mrna_data<<-mrna_data[,ind_num_mrna]
      clean_prot_data<<-prot_data[,ind_prot_num]
      if (grep("\\.[0-3]",colnames(prot_data),perl = T)){
          days<-gsub("\\.[0-3]","",colnames(prot_data)[grep("[0-9]+",colnames(prot_data))])
      }
      else{
        days<-gsub("\\.[0-3]","",colnames(mrna_data)[grep("[0-9]+",colnames(mrna_data))])
      }
      days_kiwi<<-as.numeric(as.character(days))
      total_data<-merge(mrna_data,prot_data)
      lista<-vector("list",nrow(mrna_data))
      for (i in seq(1,nrow(total_data))){
        lista[[i]]<-list("Protein_ID"=total_data[i,"Protein"],"Transcrit_ID"=total_data[i,"Transcrit"],"Transcrit_val"=as.matrix(total_data[i,3:29]),"Protein_val"=as.matrix(total_data[i,30:ncol(total_data)]),"DPA"=days_kiwi)
      }
      # test_list<<-lista
      test_list<<-sample(lista,100)
      test_el<<-sample(test_list,1)[[1]]
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
    parList<-callModule(paramList,"params",input$method_we)
    parList<<-parList$parms
    boundList<<-parList$bounds})
  observeEvent(input$fit_op,{
    print(parList())
    print(boundList())
    inFile<-input$weight_data
    poids_data<-loadData(inFile$datapath,"","",poids=T)
    print("Fitting...")
    tryCatch({
      coefs_poids<<-fitPoids_v2(poids_data[,1],poids_data[,2],input$method_we,parList(),boundList())
    },
    warning = function(warn){
      showNotification(paste0(warn), type = 'warning')
    },
    error = function(err){
      showNotification(paste0(err), type = 'err')
    })
    val_mu<-mu(c(poids_data$DPA),input$method_we,coefs_poids$coefs,coefs_poids$formula,dpa_analyse = NULL)
    data_mu<-data.frame("DPA"=c(poids_data$DPA),"Mu"=val_mu)
    g_mu<<-ggplot(data_mu,aes(x=DPA,y=Mu))+geom_line()+theme+xlab("DPA")+ylab(bquote("Growth rate "~(days^-1)))
    data_rel_mu<-data.frame("DPA"=c(poids_data$DPA),"RGR"=val_mu/fitted(coefs_poids$formula))
    g_rel_mu<<-ggplot(data_rel_mu,aes(x=DPA,y=RGR))+geom_line()+theme+xlab("DPA")+ylab(bquote("Relative growth rate "~(days^-1)))
    fit_op$state<-TRUE
    print("Finished!!")
    output$fitplot<-renderPlot({
      req((fit_op$state)==TRUE,exists("coefs_poids"))
      coefs_poids$graph
    })
    output$rates<-renderPlot({
      req((fit_op$state)==TRUE,exists(c("g_mu","g_rel_mu")))
      ggarrange(g_mu,g_rel_mu,ncol = 2)   
    })
    output$errWe<-renderText({
      paste("Error value: ",coefs_poids$error,sep="")
    })
  })
  observeEvent(input$clearGraph,{
    output$fitplot<-renderPlot(NULL)
    output$rates<-renderPlot(NULL)
    output$errWe<-renderText(NULL)
    rm("parList","boundList","g_mu","g_rel_mu","coefs_poids",envir = .GlobalEnv)
  })
  observeEvent(input$testMRNA,{
    req(exists("test_el"))
    fitR<-input$testFit_mrna
    norm_data<-normaMean(test_el$Protein_val,test_el$Transcrit_val,10)
    test_fitmrna<<-fit_testRNA(test_el$DPA,norm_data$mrna,fitR)
    output$testmrnaplot<-renderPlot({plotFitmRNA(test_el$DPA,norm_data$mrna,solmRNA(test_el$DPA,test_fitmrna$coefs,fitR))})
    output$errMrna<-renderText({
      mes1<-paste("Chosen pair: ",test_el[["Transcrit_ID"]],sep = "")
      mes2<-paste("Error value: ",test_fitmrna$error,sep = "")
      paste('<B>Error</B>')
      paste(mes1,mes2,sep = "\n")
    })
  })
  observe({
    if (input$mrnaStep=="Main calculation"){
      output$testmrnaplot<-NULL
      output$errMrna<-NULL
    }
  })
  observeEvent(input$run_loop,{
    if (input$fit_mrna!=""){
      ksmin<<-as.numeric(as.character(input$ksmin))
      score=0
      cont<-0
      poids_coef<<-coefs_poids$coefs
      formula_poids<<-coefs_poids$formula
      fitR<<-input$fit_mrna
      fitWe<<-input$method_we
      cl1 <- makeCluster(detectCores() - 1)
      # registerDoParallel(cl1)
      mess<-showNotification(paste("Running..."),duration = NULL,type = "message")
      clusterEvalQ(cl1, {
        ## set up each worker.  Could also use clusterExport()
        source("global.r")
        library(shiny)
        library(ggplot2)
        library(grid)
        library(egg)
        library(car)
        NULL
      })
      
      clusterExport(cl1,c("coefs_poids","formula_poids","ksmin","fitR","fitWe","poids_coef"))
      run_calc$run<-TRUE
      if(Sys.info()["sysname"]=="Windows"){
        num_cor<-cl1
      }
      else if(Sys.info()["sysname"]=="Linux"){
        num_cor<-detectCores()-1
      }
     # res_list<-for (el in test_list){
    res_list<-pblapply(X=test_list,function(el){
        # tryCatch({
          norm_data<-normaMean(el$Protein_val,el$Transcrit_val,ksmin)
          fittedmrna<<-fit_testRNA(el$DPA,norm_data$mrna,fitR)
          el$errorMrna<-fittedmrna$error
          el$plot_mrna<-plotFitmRNA(el$DPA,norm_data$mrna,solmRNA(el$DPA,fittedmrna$coefs,fitR))
          el$errorWeight<-coefs_poids$error
          par_k<-solgss_Borne(el$DPA,as.vector(norm_data$prot),as.numeric(norm_data$ks),score)
          if (!is.null(par_k)){
            par_k[["plot_fit_prot"]]<-plotFitProt(el$DPA,as.vector(norm_data$prot),par_k$prot_fit)
            X<-matrice_sens(el$DPA,par_k[["solK"]][,1])
            diff<-(par_k[["error"]][["errg"]][1]*norm(as.vector(norm_data$prot),"2"))^2
            par_k[["corr_matrix"]]<-matrice_corr(X,length(norm_data$prot),diff)
            el$SOL<-par_k
            # para_min<-fminunc(par_k[["solK"]][,1],fn=minSquares,time=el$DPA,exp_data=as.vector(norm_data$prot))
            
            el$confEllipse<-confidenceEllipse(el[["SOL"]][["modelList"]][["model1"]],which.coef = c("ks","kd"),fill = T,segments = 50,levels=c(0.9,0.75,0.5))
            # el$confEllipse_2<-ellipse::ellipse(el[["SOL"]][["modelList"]][["model1"]],which = c("ks","kd"),levels=c(0.9,0.75,0.5))        
            if (any(el$confEllipse[["0.9"]]<0 | el$confEllipse[["0.75"]]<0 | el$confEllipse[["0.5"]]<0)){
              el$validEllipse<-FALSE
              el$confEllipsePlot<-ggplot()+geom_polygon(as.data.frame(el[["confEllipse"]][["0.9"]]),mapping=aes(x,y),colour="red",alpha=0.3,fill="red")+geom_polygon(as.data.frame(el[["confEllipse"]][["0.75"]]),mapping=aes(x,y),colour="green",alpha=0.3,fill="green")+geom_polygon(as.data.frame(el[["confEllipse"]][["0.5"]]),mapping=aes(x,y),colour="blue",alpha=0.3,fill="blue")+theme+coord_cartesian(xlim=(c(0,max(as.data.frame(el[["confEllipse"]][["0.9"]])$x))),ylim = (c(0,max(as.data.frame(el[["confEllipse"]][["0.9"]])$y))))+ylab("kd")+xlab("ks")
            } 
            else {
              el$validEllipse<-TRUE
              el$confEllipsePlot<-ggplot()+geom_polygon(as.data.frame(el[["confEllipse"]][["0.9"]]),mapping=aes(x,y),colour="red",alpha=0.3,fill="red")+geom_polygon(as.data.frame(el[["confEllipse"]][["0.75"]]),mapping=aes(x,y),colour="green",alpha=0.3,fill="green")+geom_polygon(as.data.frame(el[["confEllipse"]][["0.5"]]),mapping=aes(x,y),colour="blue",alpha=0.3,fill="blue")+theme+ylab("kd")+xlab("ks")
              ## check ellipse package 
              }
          }
          return(el)
        # },error=function(e){showNotification(paste0("Protein fitting not achieved for ",el$Transcrit_ID,sep=" "),type = "error",duration = NULL)})
      },cl=num_cor)
      stopCluster(cl1)
      valid_res<<-Filter(function(x) {length(x) > 6}, res_list)
      mess<-showNotification(paste("Finished!!"),duration = NULL,type = "message")
      en_but$enable<-TRUE
      names_trans<-sapply(valid_res,with,Transcrit_ID)
      output$results<-renderUI({
        tagList(sidebarLayout(
          sidebarPanel(
            downloadButton("downFile","Download data table"),
            selectInput("res_trans","Select mRNA ID",choices = names_trans),
            textInput("err_th","Select error threshold",value=0.3),
            downloadButton("validTable","Download valid table")
          ),
          mainPanel(resultsKsKdUI("res_trans"))
        ))
      })
      
      all_res<-lapply(valid_res, function(el){
        res<-list()
        res[["TranscritID"]]<-el[["Transcrit_ID"]]
        res[["Weight formula"]]<-input$method_we
        res[["Weight error"]]<-el[["errorWeight"]]
        res[["mRNA formula"]]<-fitR
        res[["mRNA error"]]<-el[["errorMrna"]]
        res[["Mean mRNA concentration"]]<-mean(el[["Transcrit_val"]],na.rm = T)
        res[["Mean protein concentration"]]<-mean(el[["Protein_val"]],na.rm = T)
        res[["Starting protein concentration value"]]<-unname(el[["SOL"]][["solK"]][1,1])
        res[["Normalized ks"]]<-unname(el[["SOL"]][["solK"]][2,1])
        res[["ks"]]<-unname(el[["SOL"]][["solK"]][2,1])*res[["Mean protein concentration"]]/res[["Mean mRNA concentration"]]
        res[["kd"]]<-unname(el[["SOL"]][["solK"]][3,1])
        res[["Fitting error value"]]<-el[["SOL"]][["error"]][["errg"]][1]
        res[["Fitting error score"]]<-el[["SOL"]][["error"]][["score"]]
        res[["Fitting error message"]]<-el[["SOL"]][["error"]][["message"]]
        res[["Optimization error score"]]<-el[["SOL"]][["opt_eval"]][["score"]]
        res[["Optimization error message"]]<-el[["SOL"]][["opt_eval"]][["message"]]
        res[["Valid confidence ellipse"]]<-el[["validEllipse"]]
        return(res)
      })
      final_table<<-rbindlist(all_res)
      output$downFile<-downloadHandler(
        filename = function(){
          paste("results_KsKd-",Sys.Date(),".csv",sep="")
        },
        content = function(file){
          write.csv(final_table,file)
        },contentType = "text/csv"
        
      )
      
      
    }
  })
  
  observe({
    
    req(exists("final_table"))
    # if (as.numeric(as.character(input$err_th))<=0){
    #   showNotification("Invalid threshold value",type = "error")
    # }
    valid_table<<-final_table[which((final_table$`Fitting error value`<as.numeric(as.character(input$err_th))) & (final_table$`Valid confidence ellipse`==TRUE)),]
    output$validTable<-downloadHandler(
      filename = function(){
        paste("validResults_KsKd-",Sys.Date(),".csv",sep="")
      },
      content = function(file){
        write.csv(valid_table,file)
      },contentType = "text/csv"
      
    )
    if (nrow(valid_table)==0){
      shinyjs::disable("validTable")
      showNotification("No valid results are found",type = "error")
    }
  })
# (final_table$`Optimization error score`==10) &
  observeEvent(input$res_trans,{
    callModule(resultsKsKd,"res_trans",valid_res,input$res_trans)})
  onStop(function() {
    rm(list=ls())
  })
}