rm(list = ls())

# Create Shiny object

library(shiny)
library(shinythemes)

formula_tabs<-tabsetPanel(
  tabPanel("double_sig",
           withMathJax("$$y=d+\\frac{a}{1+exp^{-b*(t-c)}}+\\frac{e}{1+exp^{-f*(t-g)}}$$")
  ),
  tabPanel("gompertz",
           withMathJax("$$y=b*exp^{\\ln(\\frac{c}{b})*exp^{-a*t}}$$")
  ),
  tabPanel("verhulst",
           withMathJax("$$y=\\frac{b*c}{b+(b-c)*exp^{-a*t}}$$")
           
  ),
  tabPanel("log_poly",
           withMathJax("")
  ),
  
  id = "formulas",
  type = "tabs"

)

fluidPage(theme = shinytheme("united"),tags$style("#params { display:none; } #formulas { display:none; }"),
          navbarPage("Protein turnover model",id="main",
                     tabPanel("Main",
                     tabsetPanel(id="tabs",tabPanel("Input data",
                                                    sidebarLayout(
                                                      sidebarPanel(
                                                        checkboxInput("multiple","Single file",value = FALSE),
                                                        uiOutput("mult_files"),
                                                        uiOutput("sing_file"),
                                                        actionButton("disp_distr","Show distributions")
                                                      ),
                                                      mainPanel(
                                                        plotOutput("distr_plot"),
                                                        plotOutput("distr_stade")
                                                      )
                                                    )
                     ),
                     tabPanel("Weight fitting",
                              sidebarLayout(
                                sidebarPanel(
                                  tags$head(tags$style(HTML('
                                                        .selectize-input {white-space: nowrap}
                                                        #method_we+ div>.selectize-input{width: 5cm !important;}
                                                                                '))),
                                  fileInput("weight_data","Choose weight data to be fitted",accept = c("text/csv")),
                                  div(style="display:inline-block",selectInput("method_we","Select fitting formula",choices = c("Logistic"="verhulst","Gompertz"="gompertz","Empiric"="empirique","Log polynomial"="log_poly","Double sigmoid"="double_sig"))),
                                  div(style="display:inline-block",formula_tabs), ## "Contois"="contois",,"Noyau"="seed",
                                  paramListInput("params"),
                                  div(style="display:inline-block",actionButton("fit_op","Fit")),
                                  div(style="display:inline-block",actionButton("clearGraph","Clear graphs"))
                                ),
                                mainPanel(
                                  plotOutput("fitplot"),
                                  plotOutput("rates"),
                                  textOutput("errWe")
                                )
                              )
                              
                     ),
                     tabPanel("mRNA fitting and calculation",
                              sidebarLayout(
                                sidebarPanel(
                                  tabsetPanel(type="pills",id="mrnaStep",
                                              tabPanel("Main calculation",textInput("ksmin","Value of ksmin",value =3*4*3*3.6*24),
                                                       selectInput("fit_mrna","Select fitting formula",choices=c("3rd degree polynomial"="3_deg","6th degree polynomial"="6_deg","3rd degree logarithmic polynomial"="3_deg_log")),
                                                       actionButton("run_loop","Run calculation")
                                                       ),
                                              tabPanel("Test mRNA fitting",
                                                       selectInput("testFit_mrna","Select fitting formula",choices=c("3rd degree polynomial"="3_deg","6th degree polynomial"="6_deg","3rd degree logarithmic polynomial"="3_deg_log")),
                                                       actionButton("testMRNA","Run mRNA fitting test"))
                                  )
                                ),
                                mainPanel(plotOutput("testmrnaplot",width = "100%"),
                                          textOutput("errMrna"))
                              )
                              
                     ),
                     tabPanel("Results",id="tabRes",
                              uiOutput("results")
                              
                              
                     ))),
                     tabPanel("Help", 
                              includeMarkdown("tutorial.md"),)
                     
                    )
)