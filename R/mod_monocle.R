mod_Monocle_ui <- function(id) {
  ns <- NS(id)
  
  # Monocle
  tabItem(tabName = "monocle",
          navbarPage("Analyses: ",
                     tabPanel("Trajectory", 
                              radioGroupButtons(inputId = "color_trajectory", label = "Color by :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                              fluidRow(),
                              plotOutput("cell_trajectory", width = "100%",  height = "600px")
                     ),
                     tabPanel("Ordering", 
                              radioGroupButtons(inputId = "color_order", label = "Color by :", choices = singlet@tools$meta_variable, justified = TRUE, checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                              fluidRow(),
                              selectInput('DDvariables', 'Choices', NULL, selectize=TRUE, multiple=TRUE, selected = NULL),
                              uiOutput("DDorder_trajectory")
                     )
                     
          ),
  )
}

mod_Monocle_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_Monocle_ui("Monocle_ui_1")

## To be copied in the server
# callModule(mod_Monocle_server, "Monocle_ui_1")