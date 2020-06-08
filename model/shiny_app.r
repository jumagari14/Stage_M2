# Load packages
# # Load data
# trend_data <- read_xlsx("Paires_mrna_prot_kiwi_nouvMW.xlsx",sheet = "Proteines")
# trend_description <- read_xlsx("Paires_mrna_prot_kiwi_nouvMW.xlsx",sheet = "Transcrits")
# 
# # Define UI
# ui <- fluidPage(theme = shinytheme("lumen"),
#                 titlePanel("Protein turnover model"),
#                 sidebarLayout(
#                   sidebarPanel(
#                     
#                     # Select type of trend to plot
#                     selectInput(inputId = "type", label = strong("Trend index"),
#                                 choices = unique(trend_data$type),
#                                 selected = "Travel"),
#                     
#                     # Select date range to be plotted
#                     dateRangeInput("date", strong("Date range"), start = "2007-01-01", end = "2017-07-31",
#                                    min = "2007-01-01", max = "2017-07-31"),
#                     
#                     # Select whether to overlay smooth trend line
#                     checkboxInput(inputId = "smoother", label = strong("Overlay smooth trend line"), value = FALSE),
#                     
#                     # Display only if the smoother is checked
#                     conditionalPanel(condition = "input.smoother == true",
#                                      sliderInput(inputId = "f", label = "Smoother span:",
#                                                  min = 0.01, max = 1, value = 0.67, step = 0.01,
#                                                  animate = animationOptions(interval = 100)),
#                                      HTML("Higher values give more smoothness.")
#                     )
#                   ),
#                   
#                   # Output: Description, lineplot, and reference
#                   mainPanel(
#                     plotOutput(outputId = "lineplot", height = "300px"),
#                     textOutput(outputId = "desc"),
#                     tags$a(href = "https://www.google.com/finance/domestic_trends", "Source: Google Domestic Trends", target = "_blank")
#                   )
#                 )
# )
# 
# # Define server function
# server <- function(input, output) {
#   
#   # Subset data
#   selected_trends <- reactive({
#     req(input$date)
#     validate(need(!is.na(input$date[1]) & !is.na(input$date[2]), "Error: Please provide both a start and an end date."))
#     validate(need(input$date[1] < input$date[2], "Error: Start date should be earlier than end date."))
#     trend_data %>%
#       filter(
#         type == input$type,
#         date > as.POSIXct(input$date[1]) & date < as.POSIXct(input$date[2]
#         ))
#   })
#   
#   
#   # Create scatterplot object the plotOutput function is expecting
#   output$lineplot <- renderPlot({
#     color = "#434343"
#     par(mar = c(4, 4, 1, 1))
#     plot(x = selected_trends()$date, y = selected_trends()$close, type = "l",
#          xlab = "Date", ylab = "Trend index", col = color, fg = color, col.lab = color, col.axis = color)
#     # Display only if smoother is checked
#     if(input$smoother){
#       smooth_curve <- lowess(x = as.numeric(selected_trends()$date), y = selected_trends()$close, f = input$f)
#       lines(smooth_curve, col = "#E6553A", lwd = 3)
#     }
#   })
#   
#   # Pull in description of trend
#   output$desc <- renderText({
#     trend_text <- filter(trend_description, type == input$type) %>% pull(text)
#     paste(trend_text, "The index is set to 1.0 on January 1, 2004 and is calculated only for US search traffic.")
#   })
# }

# Create Shiny object

library(shiny)
library(shinythemes)
library(dplyr)
library(readr)
library(readxl)
library(shinyjs)
library("gplots")
library("reshape2")
library("vegan")
library("ggplot2")
library("grid")
library("goeveg")
library(gridExtra)
library(cowplot)
library(egg)

# source("input.r")
source("functions.r")

ui <- fluidPage(useShinyjs(),theme = shinytheme("lumen"),
                # singleton(tags$head(HTML(
                #   '
                #                         <script type="text/javascript">
                #                         $(document).ready(function() {
                #                         // disable download at startup. data_file is the id of the downloadButton
                #                         $("#data_file").attr("disabled", "true").attr("onclick", "return false;");
                #                         
                #                         Shiny.addCustomMessageHandler("download_ready", function(message) {
                #                         $("#data_file").removeAttr("disabled").removeAttr("onclick").html(
                #                         "<i class=\\"fa fa-download\\"></i>Download (file size: " + message.fileSize + ")");
                #                         });
                #                         })
                #                         </script>
                #                         '
                # ))),
                navbarPage("Protein turnover model",id="tabs",
                           tabsetPanel(tabPanel("Input data",
                                    checkboxInput("multiple","Single file",value = FALSE),
                                    uiOutput("mult_files"),
                                    uiOutput("sing_file"),
                                    actionButton("disp_distr","Show distributions"),
                                    plotOutput("distr_plot"),
                                    plotOutput("distr_stade")
                            ),
                           tabPanel("Weight fitting",
                                    fileInput("weight_data","Choose weight data to be fitted",accept = c("text/csv")),
                                    selectInput("method_we","Select fitting formula",choices = c("Logistic"="verhulst","Gompertz"="gompertz","Contois"="contois","Empiric"="empirique","Noyau"="seed","Log polynomial"="log_poly","Double sigmoid"="double_sig")),
                                    uiOutput("doubl_sig"),
                                    uiOutput("gompertz"),
                                    actionButton("fit_op","Fit"),
                                    plotOutput("fitplot")
                            ),
                           tabPanel("mRNA fitting and calculation",
                                    textInput("ksmin","Value of ksmin",value =3*4*3*3.6*24),
                                    selectInput("fit_mrna","Select fitting formula",choices=c("3rd degree polynomial"="3_deg","6th degree polynomial"="6_deg","3rd degree logarithmic polynomial"="3_deg_log")),
                                    actionButton("run_loop","Run calculation"),
                                    disabled(downloadButton("downFile","Save results")))
                                       ))
  # useShinyjs()
  # uiOutput("header"),
  # br(),
  # actionButton("prevBtn", "< Previous"),
  # actionButton("nextBtn", "Next >")
                )


server <- function(input, output, session) {
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
    
  output$doubl_sig<-renderUI({
    if (input$method_we=="double_sig"){
    # print("Jaja")
    tagList(textInput("par1_sig","Enter value of par1",value = 48),
            textInput("par2_sig","Enter value of par2",value = 0.144),
            textInput("par3_sig","Enter value of par3",value = 35),
            textInput("par4_sig","Enter value of par4",value = 0.4),
            textInput("par5_sig","Enter value of par5",value = 48),
            textInput("par6_sig","Enter value of par6",value = 0.042),
            textInput("par7_sig","Enter value of par7",value = 90))
    }
  })
  output$gompertz<-renderUI({
    if (input$method_we=="gompertz"){
        tagList(textInput("par1_sig","Enter value of par1",value =0.065),
                textInput("par2_sig","Enter value of par2",value = 114.39),
                textInput("par3_sig","Enter value of par3",value = 0.52))
      }
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
      x_ind<-grep("par*",names(x))
      newlist<-vector("list",length(x_ind))
      names(newlist)<-names(x[x_ind])
      for (el in names(newlist)){
        newlist[[el]]<-as.numeric(as.character(input[[el]]))
      }
      names(newlist)<-gsub("_sig","",names(newlist))
      newlist<-newlist[order(names(newlist))]
      newlist
    })
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
          if (!is.null(par_k)){
            test_list[[cont]]$SOL<-par_k
            print(test_list[[cont]]$SOL)
            # write.csv(test_list[[cont]][["SOL"]][["solK"]],paste("solK/",paste(test_list[[cont]][["Transcrit_ID"]],"_Sol_ks_kd.csv"),sep = ""))
          }
        },error=function(e){showNotification(paste0("Protein fitting not achieved for ",el$Transcrit_ID,sep=" "),type = "error",duration = NULL)})
        
      }
      valid_res<<-Filter(function(x) {length(x) > 5}, test_list)
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
      print("Prueba")
      enable("downFile")
    }
  })
  }


if (interactive()) shinyApp(ui = ui, server = server)

