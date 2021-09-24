mod_Metadata_ui <- function(id) {
  ns <- NS(id)
  # Métadonnées
  tabItem(tabName = "metadata",
          h1("Présentation du patient et métadonnées : "),
          splitLayout(cellWidths=c("35%","65%"),verticalLayout(fluid = T,
                                                               tags$head(tags$style(HTML('.shiny-split-layout>div {overflow: hidden;}')),),
                                                               verbatimTextOutput("nb_clonotype_cell2"),br(),br(),
                                                               HTML('<center><img src="bcr.png" width="400"></center>')),
                      plotlyOutput('dataTable', width = "100%",  height = "800px")
          ),
          
          downloadButton("PDF_report", "PDF report"),
          downloadButton("HTML_report", "HTML report"),
          downloadButton("Word_report", "Word report")
  )}

mod_Metadata_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_Metadata_ui("Metadata_ui_1")

## To be copied in the server
# callModule(mod_Metadata_server, "Metadata_ui_1")