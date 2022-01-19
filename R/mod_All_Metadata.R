mod_All_Metadata_ui <- function(id) {
  ns <- NS(id)
  # metadata Metadata
  tabItem(tabName = "All_Metadata",
    navbarPage("Information : ",
      tabPanel("Phenotype",
        fluidRow(
          box("", width = 4, plotlyOutput(ns("FL08_Meta"), width = "100%",  height = "600px")),
          box("", width = 4, plotlyOutput(ns("FL14_Meta"), width = "100%",  height = "600px")),
          box("", width = 4, plotlyOutput(ns("FL12_Meta"), width = "100%",  height = "600px"))),
        fluidRow(
          box("", width = 4, plotlyOutput(ns("FL09_Meta"), width = "100%",  height = "600px")),
          box("", width = 4, plotlyOutput(ns("FL02_Meta"), width = "100%",  height = "600px")),
          box("", width = 4, plotlyOutput(ns("FL05_Meta"), width = "100%",  height = "600px"))),          
        fluidRow(
          box("", width = 4, plotlyOutput(ns("FL082_Meta"), width = "100%",  height = "600px")),
          box("", width = 4, plotlyOutput(ns("FL06_Meta"), width = "100%",  height = "600px")))
      ),
      tabPanel("VDJ",
          fluidRow(
            box("", width = 4, plotlyOutput(ns("FL08_VDJ"), width = "100%",  height = "600px")),
            box("", width = 4, plotlyOutput(ns("FL14_VDJ"), width = "100%",  height = "600px")),
            box("", width = 4, plotlyOutput(ns("FL12_VDJ"), width = "100%",  height = "600px"))),
          fluidRow(
            box("", width = 4, plotlyOutput(ns("FL02_VDJ"), width = "100%",  height = "600px")),
            box("", width = 4, plotlyOutput(ns("FL05_VDJ"), width = "100%",  height = "600px")),
            box("", width = 4, plotlyOutput(ns("FL09_VDJ"), width = "100%",  height = "600px"))),
        fluidRow(
          box("", width = 4, plotlyOutput(ns("FL082_VDJ"), width = "100%",  height = "600px")),
          box("", width = 4, plotlyOutput(ns("FL06_VDJ"), width = "100%",  height = "600px"))))
      )
  )
}

mod_All_Metadata_server <- function(input, output, session) {
  ns <- session$ns
  
  ## -- Méta-Phenotpe -- ##
  output$FL08_Meta = renderPlotly({load(file = "inst/app/www/FL08G0293/FL08G0293_VDJ.RData") ; f})
  output$FL09_Meta = renderPlotly({load(file = "inst/app/www/FL09C1164/FL09C1164_VDJ.RData") ; f})
  output$FL12_Meta = renderPlotly({load(file = "inst/app/www/FL12C1888/FL12C1888_VDJ.RData") ; f})
  output$FL14_Meta = renderPlotly({load(file = "inst/app/www/FL140304/FL140304_VDJ.RData") ; f})
  output$FL02_Meta = renderPlotly({load(file = "inst/app/www/FL02G095/FL02G095_VDJ.RData") ; f})
  output$FL05_Meta = renderPlotly({load(file = "inst/app/www/FL05G0330/FL05G0330_VDJ.RData") ; f})
  output$FL082_Meta = renderPlotly({load(file = "inst/app/www/FL08G0431/FL08G0431_VDJ.RData") ; f})
  output$FL06_Meta = renderPlotly({load(file = "inst/app/www/FL06G1206/FL06G1206_VDJ.RData") ; f})
  
  ## -- Méta-VDJ -- ##
  output$FL08_VDJ = renderPlotly({load(file = "inst/app/www/FL08G0293/FL08G0293_VDJ.RData") ; e})
  output$FL09_VDJ = renderPlotly({load(file = "inst/app/www/FL09C1164/FL09C1164_VDJ.RData") ; e})
  output$FL12_VDJ = renderPlotly({load(file = "inst/app/www/FL12C1888/FL12C1888_VDJ.RData") ; e})
  output$FL14_VDJ = renderPlotly({load(file = "inst/app/www/FL140304/FL140304_VDJ.RData") ; e})
  output$FL02_VDJ = renderPlotly({load(file = "inst/app/www/FL02G095/FL02G095_VDJ.RData") ; e})
  output$FL05_VDJ = renderPlotly({load(file = "inst/app/www/FL05G0330/FL05G0330_VDJ.RData") ; e})
  output$FL082_VDJ = renderPlotly({load(file = "inst/app/www/FL08G0431/FL08G0431_VDJ.RData") ; e})
  output$FL06_VDJ = renderPlotly({load(file = "inst/app/www/FL06G1206/FL06G1206_VDJ.RData") ; e})
}


## To be copied in the UI
# mod_All_Metadata_ui("meta_Metadata_ui_1")

## To be copied in the server
# callModule(mod_All_Metadata_server, "meta_Metadata_ui_1")