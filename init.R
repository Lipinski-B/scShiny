# init.R
#
# Example R code to install packages if not already installed
#

my_packages = c('plotly','shiny','shinybusy','shinydashboard','dashboardthemes','shinyWidgets')#'Seurat',


install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages(p)
    BiocManager::install("dittoSeq")
  }
}

invisible(sapply(my_packages, install_if_missing))

