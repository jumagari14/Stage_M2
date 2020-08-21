list.of.packages<-"shiny"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos="https://pbil.univ-lyon1.fr/CRAN/")
lapply(list.of.packages,require,character.only=TRUE)

x <- system("ipconfig", intern=TRUE)
z <- x[grep("IPv4", x)]
ip <- gsub(".*? ([[:digit:]])", "\\1", z)
print(paste0("the Shiny Web application runs on: http://", ip))
if (length(ip)>1) ip=ip[1]
shiny::runApp('./model',host=ip[1],launch.browser = T)