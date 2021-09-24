mod_PCA_ui <- function(id) {
  ns <- NS(id)
  
  # 10 PCA dimensions
  tabItem(tabName = "PCA_facteurs",
          navbarPage("PCA paramètres : ",
                     tabPanel("Facteurs", verbatimTextOutput("PCA10")),
                     tabPanel("ElbowPlot", plotOutput("ElbowPlot", width = "100%",  height = "1000px")),
                     tabPanel("DimHeatmap", plotOutput("DimHeatmap", width = "100%",  height = "1000px")),
                     tabPanel("VizDimLoadings", plotOutput("VizDimLoadings", width = "100%",  height = "1000px")),
                     tabPanel("JackStrawPlot", plotOutput("JackStrawPlot", width = "100%",  height = "1000px"))
          ),
  )
  
  
}

mod_PCA_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_PCA_ui("PCA_ui_1")

## To be copied in the server
# callModule(mod_PCA_server, "PCA_ui_1")