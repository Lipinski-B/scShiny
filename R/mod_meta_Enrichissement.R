mod_meta_Enrichissement_ui <- function(id) {
  ns <- NS(id)
  # metadata enrichissement
  tabItem(tabName = "metadata_enrichissement",
          navbarPage("Figure : ",
                     tabPanel("Histograme",plotlyOutput("hist_enrichissement", width = "100%",  height = "1000px")),
                     tabPanel("Datatable",DT::DTOutput("df_enrichissement")))
  )
}

mod_meta_Enrichissement_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_meta_Enrichissement_ui("meta_Enrichissement_ui_1")

## To be copied in the server
# callModule(mod_meta_Enrichissement_server, "meta_Enrichissement_ui_1")