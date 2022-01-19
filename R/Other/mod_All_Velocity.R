mod_All_Velocity_ui <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "All_Velocity",
          navbarPage("FL08G0431 : ",
                     tabPanel("Avec",
                              splitLayout(cellWidths = c("50%", "50%"),
                                          HTML('<center><img src="Velocity/scVelo_all_sans_FL082.png" width="1000"></center>'),
                                          HTML('<center><img src="Velocity/scVelo_B_RCHOP_sans_FL082.png" width="1000"></center>'))),
                     tabPanel("Sans", 
                              splitLayout(cellWidths = c("50%", "50%"),                                          
                                          HTML('<center><img src="Velocity/scVelo_all.png" width="1000"></center>'),
                                          HTML('<center><img src="Velocity/scVelo_B_RCHOP.png" width="1000"></center>')))   

          )
  )
}

mod_All_Velocity_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_All_Velocity_ui("All_Velocity_ui_1")

## To be copied in the server
# callModule(mod_All_Velocity_server, "All_Velocity_ui_1")