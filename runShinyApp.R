# 
# R Script to load GUI 
# Garcia, Juan Manuel
# 25/08/2020
# 

# Load shiny package, install if necessary
list.of.packages<-c("shiny","shinythemes","shinyjs","rlang")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos="https://pbil.univ-lyon1.fr/CRAN/")
lapply(list.of.packages,require,character.only=TRUE)
# Open GUI in default web browser
shiny::runApp('./model',launch.browser = T)