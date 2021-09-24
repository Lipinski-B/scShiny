mod_VDJ_ui <- function(id) {
  ns <- NS(id)
  # VDJ
  tabItem(tabName = "VDJ",
          box("Clonotype", width = 12, plotlyOutput("VDJ_Clonotype", width = "100%",  height = "400px")),
          box("Heavy", width = 7, plotlyOutput("VDJ_Heavy", width = "100%",  height = "400px"),
              column(1,align="right",checkboxGroupInput("heavy_details", NULL, choices = list("Details" = "Details"), selected = 0)),
              conditionalPanel(
                condition = "input.heavy_details == 'Details' ", br(),
                plotlyOutput("VDJ_DHeavy", width = "100%",  height = "400px"))
          ),
          box("Lights", width = 5, plotlyOutput("VDJ_Light", width = "100%",  height = "400px"),
              column(1,align="right",checkboxGroupInput("light_details", NULL, choices = list("Details" = "Details"), selected = 0)),
              conditionalPanel(
                condition = "input.light_details == 'Details' ", br(),
                plotlyOutput("VDJ_DLight", width = "100%",  height = "400px"))
          ),
          fluidRow(
            box("V", width = 4, plotlyOutput("V", width = "100%",  height = "400px")),
            box("D", width = 4, plotlyOutput("D", width = "100%",  height = "400px")),
            box("J", width = 4, plotlyOutput("J", width = "100%",  height = "400px")))
  )
}

mod_VDJ_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_VDJ_ui("VDJ_ui_1")

## To be copied in the server
# callModule(mod_VDJ_server, "VDJ_ui_1")