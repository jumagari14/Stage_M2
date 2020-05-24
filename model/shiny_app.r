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

source("main.r")
source("input.r")
source("../../Analyse_stats/functions.r")

ui <- fluidPage(theme = shinytheme("lumen"),
                navbarPage("Protein turnover model",id="tabs",
                           tabsetPanel(tabPanel("Input data",
                                    checkboxInput("multiple","Single file",value = FALSE),
                                    uiOutput("mult_files"),
                                    uiOutput("sing_file"),
                                    actionButton("disp_distr","Show distributions"),
                                    plotOutput("distr_plot")
                                    ),
                           tabPanel("Weight fitting",
                                    fileInput("weight_data","Choose weight data to be fitted",accept = c("text/csv")),
                                    selectInput("method_we","Select fitting formula",choices = c("Logistic"="logis","Gompertz"="gomp","Contois"="conto","Empiric"="empi","Noyau"="seed","Log polynomial"="logply","Double sigmoid"="doubl_sig")),
                                    actionButton("fit_op","Fit"),
                           ))
  # useShinyjs()
  # uiOutput("header"),
  # br(),
  # actionButton("prevBtn", "< Previous"),
  # actionButton("nextBtn", "Next >")
                )
)

server <- function(input, output, session) {
  # rv <- reactiveValues(page = 1)
  # output$header<-renderText({ifelse(rv$page==2,"Input data","Second step")})
  # observe({
  #   toggleState(id = "prevBtn", condition = rv$page > 1)
  #   toggleState(id = "nextBtn", condition = rv$page < NUM_PAGES)
  #   hide(selector = ".page")
  #   show(paste0("step", rv$page))
  # })

  # navPage <- function(direction) {
  #   rv$page <- rv$page + direction
  # }

  # observeEvent(input$prevBtn, navPage(-1))
  # observeEvent(input$nextBtn, navPage(1))
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
  observe({
    if(!is.null(input$data_file)){
      inFile<-input$data_file
      list_data<-loadData(inFile$datapath,input$rna_tab,input$protein_tab)
      mrna_data<-list_data$mrna
      prot_data<-list_data$prot
      clean_mrna_data<-mrna_data[,-which(is.na(as.numeric(as.character(colnames(mrna_data)))))]
      clean_prot_data<-prot_data[,-which(is.na(as.numeric(as.character(colnames(prot_data)))))]
      
      observeEvent(input$disp_distr,{
        print("Plotting...")
        output$distr_plot<-renderPlot({print(combineGraphs(clean_mrna_data,clean_prot_data,"",moyenne = T))})
        print("Finished")
        })

    }
  })
  
}

if (interactive()) shinyApp(ui = ui, server = server)

