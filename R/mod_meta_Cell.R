mod_meta_Cell_ui <- function(id) {
  ns <- NS(id)
  # metadata cell
  tabItem(tabName = "metadata_cell",
          fluidRow(
            box("", width = 6, plotlyOutput("FL08_Meta", width = "100%",  height = "600px")),
            box("", width = 6, plotlyOutput("FL09_Meta", width = "100%",  height = "600px"))),
          fluidRow(
            box("", width = 6, plotlyOutput("FL12_Meta", width = "100%",  height = "600px")),
            box("", width = 6, plotlyOutput("FL14_Meta", width = "100%",  height = "600px"))),
          fluidRow(
            box("", width = 6, plotlyOutput("FL02_Meta", width = "100%",  height = "600px")),
            box("", width = 6, plotlyOutput("FL05_Meta", width = "100%",  height = "600px"))),
  )
}

mod_meta_Cell_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_meta_Cell_ui("meta_Cell_ui_1")

## To be copied in the server
# callModule(mod_meta_Cell_server, "meta_Cell_ui_1")