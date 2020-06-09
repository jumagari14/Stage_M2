library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(grid)
library(egg)

# source("input.r")
source("functions.r")

function(input, output, session) {
  js$disableTab("tabRes")
  fit_op<-reactiveValues(data=NULL)
  parList<-reactiveValues(data=NULL)
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
  
  # output$print_form<-renderUI({
  #   if (input$method_we=="double_sig"){
  #     browser()
  #     withMathJax("$$y=d+\\frac{a}{1+exp^{-b*(t-c)}}+\\frac{e}{1+exp^{-f*(t-g)}}$$")
  #   }
  #   if (input$method_we=="gompertz"){
  #     withMathJax("$$y=b*exp^{\\ln(\\frac{par3}{par2})*exp^{-a*t}}$$")
  #   }
  #   if (input$method_we=="verhulst"){
  #     withMathJax("$$y=\\frac{b*c}{b+(b-c)*exp^{-a*t}}$$")
  #   }
  # })
  observeEvent(input$method_we, {
    updateTabsetPanel(session, "params", selected = input$method_we)
    updateTabsetPanel(session,"formulas",selected = input$method_we)
  })
  # output$doubl_sig<-renderUI({
  #   if (input$method_we=="double_sig"){
  #   # print("Jaja")
  #   tagList(textInput("par1_sig","Enter value of a",value = 48),
  #           textInput("par2_sig","Enter value of b",value = 0.144),
  #           textInput("par3_sig","Enter value of c",value = 35),
  #           textInput("par4_sig","Enter value of d",value = 0.4),
  #           textInput("par5_sig","Enter value of e",value = 48),
  #           textInput("par6_sig","Enter value of f",value = 0.042),
  #           textInput("par7_sig","Enter value of g",value = 90))
  #   }
  # })
  # output$gompertz<-renderUI({
  #   if (input$method_we=="gompertz"){
  #       tagList(textInput("par1_sig","Enter value of a",value =0.065),
  #               textInput("par2_sig","Enter value of b",value = 114.39),
  #               textInput("par3_sig","Enter value of c",value = 0.52))
  #     }
  #   })
  
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
  
  
  observeEvent(input$fit_op,{
    parList<<-reactive({
      x<-reactiveValuesToList(input)
      x_ind<-grep("par[1-9]+",names(x),perl = T)
      newlist<-vector("list",length(x_ind))
      names(newlist)<-names(x[x_ind])
      for (el in names(newlist)){
        newlist[[el]]<-as.numeric(as.character(input[[el]]))
      }
      names(newlist)<-gsub("_sig","",names(newlist))
      newlist<-newlist[order(names(newlist))]
      newlist
    })
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