require(shiny)

x <- system("ipconfig", intern=TRUE)
z <- x[grep("IPv4", x)]
ip <- gsub(".*? ([[:digit:]])", "\\1", z)
print(paste0("the Shiny Web application runs on: http://", ip))
if (length(ip)>1) ip=ip[1]
shiny::runApp('./model',host=ip[1],launch.browser = T)