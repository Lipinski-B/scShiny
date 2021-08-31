# init.R
#
# Example R code to install packages if not already installed
#

my_packages = c("Seurat", "monocle","escape","dittoSeq","plotly","dplyr","shiny","shinyjs","shinybusy","shinyWidgets","shinydashboard","dashboardthemes","future")


install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages(p)
    BiocManager::install(p)
  }
}

invisible(sapply(my_packages, install_if_missing))