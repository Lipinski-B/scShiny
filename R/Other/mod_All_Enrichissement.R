mod_All_Enrichissement_ui <- function(id) {
  ns <- NS(id)
  # metadata enrichissement
  tabItem(tabName = "All_enrichissement",
          navbarPage("Figure : ",
                     tabPanel("Histograme",plotlyOutput(ns("hist_enrichissement"), width = "100%",  height = "1000px")),
                     tabPanel("Datatable",DT::DTOutput(ns("df_enrichissement"))))
  )
}

mod_All_Enrichissement_server <- function(input, output, session) {
  ns <- session$ns
  
  output$hist_enrichissement = renderPlotly({load(file = "www/enrichissement.Rdata") ; hist})
  output$df_enrichissement =  DT::renderDT({load(file = "www/enrichissement.Rdata") ; datatable})
  
  
}


## To be copied in the UI
# mod_All_Enrichissement_ui("meta_Enrichissement_ui_1")

## To be copied in the server
# callModule(mod_All_Enrichissement_server, "meta_Enrichissement_ui_1")