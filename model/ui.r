jscode <- "
shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"

css <- "
.nav li a.disabled {
  background-color: #aaa !important;
  color: #333 !important;
  cursor: not-allowed !important;
  border-color: #aaa !important;
}"
# Create Shiny object

library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(grid)
library(egg)
library(shinythemes)

# source("input.r")

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

fluidPage(useShinyjs(),theme = shinytheme("united"),useShinyjs(),tags$style("#params { display:none; } #formulas { display:none; }"),
          extendShinyjs(text = jscode),inlineCSS(css),
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
                                  actionButton("fit_op","Fit")
                                ),
                                mainPanel(
                                  plotOutput("fitplot"),
                                  textOutput("errWe")
                                )
                              )
                              
                     ),
                     tabPanel("mRNA fitting and calculation",
                              sidebarLayout(
                                sidebarPanel(
                                  tabsetPanel(type="pills",id="mrnaStep",
                                              tabPanel("Test mRNA fitting",
                                                       selectInput("testFit_mrna","Select fitting formula",choices=c("3rd degree polynomial"="3_deg","6th degree polynomial"="6_deg","3rd degree logarithmic polynomial"="3_deg_log")),
                                                       actionButton("testMRNA","Run mRNA fitting test")),
                                              tabPanel("Main calculation",textInput("ksmin","Value of ksmin",value =3*4*3*3.6*24),
                                                       selectInput("fit_mrna","Select fitting formula",choices=c("3rd degree polynomial"="3_deg","6th degree polynomial"="6_deg","3rd degree logarithmic polynomial"="3_deg_log")),
                                                       actionButton("run_loop","Run calculation"))
                                  )
                                ),
                                mainPanel(plotOutput("testmrnaplot"),
                                          textOutput("errMrna"))
                              )
                              
                     ),
                     tabPanel("Results",id="tabRes",
                              sidebarLayout(
                                sidebarPanel(
                                  downloadButton("downFile","Download data table"),
                                  uiOutput("results")
                                ),
                                mainPanel(resultsKsKdUI("res_trans"))
                              )
                              
                              
                     ))),
                     tabPanel("Help", 
                              includeMarkdown("tutorial.md"),)
                     
                    )
)