library(shiny)

port <- Sys.getenv('PORT')

shiny::runApp(
  appDir = setwd("/home/boris/Bureau/scShiny/shiny/"),
  host = '0.0.0.0',
  port = as.numeric(port)
)