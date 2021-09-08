# init.R
#
# Example R code to install packages if not already installed
#

my_packages = c('plotly','shiny','shinybusy','shinydashboard','dashboardthemes','shinyWidgets')#'Seurat',


install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages("htmltools_0.5.2.tar.gz", repos = NULL, type="source")
    install.packages(p)
    install.packages("dittoSeq_1.4.1.tar.gz", repos = NULL, type="source")
  }
}

invisible(sapply(my_packages, install_if_missing))

