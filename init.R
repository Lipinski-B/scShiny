# init.R
#
# Example R code to install packages if not already installed
#

my_packages = c("plotly","dplyr","shiny","shinyjs","shinybusy","shinyWidgets","shinydashboard","dashboardthemes")#"Seurat") 


install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    helpers.installPackages(p)
    #BiocManager::install("monocle","escape","dittoSeq")
  }
}

invisible(sapply(my_packages, install_if_missing))