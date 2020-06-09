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

# source("input.r")
source("functions.r")

parameter_tabs <- tabsetPanel(
  tabPanel("double_sig",
           numericInput("par1_sig","Enter value of a",value = 48),
           numericInput("par2_sig","Enter value of b",value = 0.144),
           numericInput("par3_sig","Enter value of c",value = 35),
           numericInput("par4_sig","Enter value of d",value = 0.4),
           numericInput("par5_sig","Enter value of e",value = 48),
           numericInput("par6_sig","Enter value of f",value = 0.042),
           numericInput("par7_sig","Enter value of g",value = 90)
  ),
  tabPanel("gompertz",
           numericInput("par1_sig","Enter value of a",value =0.065),
           numericInput("par2_sig","Enter value of b",value = 114.39),
           numericInput("par3_sig","Enter value of c",value = 0.52)
  ),
  tabPanel("verhulst",
           numericInput("par1_sig","Enter value of a",value =0.1),
           numericInput("par2_sig","Enter value of b",value = 100),
           numericInput("par3_sig","Enter value of c",value = 1)
           
  ),
  id = "params",
  type = "tabs"
  
  # tabPanel("exponential",
  #          numericInput("rate", "rate", value = 1, min = 0),
  # )
)

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
  id = "formulas",
  type = "tabs"
  
  # tabPanel("exponential",
  #          numericInput("rate", "rate", value = 1, min = 0),
  # )
)

fluidPage(useShinyjs(),theme = shinytheme("lumen"),useShinyjs(),tags$style("#params { display:none; } #formulas { display:none; }"),
          extendShinyjs(text = jscode),inlineCSS(css),
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
                              div(style="display:inline-block",selectInput("method_we","Select fitting formula",choices = c("Logistic"="verhulst","Gompertz"="gompertz","Contois"="contois","Empiric"="empirique","Noyau"="seed","Log polynomial"="log_poly","Double sigmoid"="double_sig"))),
                              div(style="display:inline-block",formula_tabs),
                              parameter_tabs,
                              actionButton("fit_op","Fit"),
                              plotOutput("fitplot")
                     ),
                     tabPanel("mRNA fitting and calculation",
                              textInput("ksmin","Value of ksmin",value =3*4*3*3.6*24),
                              selectInput("fit_mrna","Select fitting formula",choices=c("3rd degree polynomial"="3_deg","6th degree polynomial"="6_deg","3rd degree logarithmic polynomial"="3_deg_log")),
                              actionButton("run_loop","Run calculation"),
                              disabled(downloadButton("downFile","Save results"))
                     ),
                     tabPanel("Results",id="tabRes",
                              uiOutput("select_res"),
                              plotOutput("fit_prot_plot")
                     )
                     ))
          # useShinyjs()
          # uiOutput("header"),
          # br(),
          # actionButton("prevBtn", "< Previous"),
          # actionButton("nextBtn", "Next >")
)