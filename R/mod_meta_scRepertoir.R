mod_meta_scRepertoir_ui <- function(id) {
  ns <- NS(id)
  # metadata scRepertoir
  tabItem(tabName = "metadata_scRepertoir",
          navbarPage("Figure : ",
                     tabPanel("Homeostasis",HTML('<left><img src="scRepertoir/clonalHomeostasis.png" width="1000"></left>')),
                     tabPanel("Proportion",HTML('<left><img src="scRepertoir/clonalProportion.png" width="900"></left>')),
                     tabPanel("Contig length",HTML('<left><img src="scRepertoir/lengthContig.png" width="900"></left>')),
                     tabPanel("Contig quant",HTML('<left><img src="scRepertoir/quantContig_output.png" width="700"></left>')))
  )
}

mod_meta_scRepertoir_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_meta_scRepertoir_ui("meta_scRepertoir_ui_1")

## To be copied in the server
# callModule(mod_meta_scRepertoir_server, "meta_scRepertoir_ui_1")