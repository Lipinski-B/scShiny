mod_meta_VDJ_ui <- function(id) {
  ns <- NS(id)
  # metadata VDJ
  tabItem(tabName = "metadata_VDJ",
          fluidRow(
            box("", width = 6, plotlyOutput("FL08_VDJ", width = "100%",  height = "600px")),
            box("", width = 6, plotlyOutput("FL09_VDJ", width = "100%",  height = "600px"))),
          fluidRow(
            box("", width = 6, plotlyOutput("FL12_VDJ", width = "100%",  height = "600px")),
            box("", width = 6, plotlyOutput("FL14_VDJ", width = "100%",  height = "600px"))),
          fluidRow(
            box("", width = 6, plotlyOutput("FL02_VDJ", width = "100%",  height = "600px")),
            box("", width = 6, plotlyOutput("FL05_VDJ", width = "100%",  height = "600px"))),
          
  )
}

mod_meta_VDJ_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_meta_VDJ_ui("meta_VDJ_ui_1")

## To be copied in the server
# callModule(mod_meta_VDJ_server, "meta_VDJ_ui_1")