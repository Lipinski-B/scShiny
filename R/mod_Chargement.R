mod_Chargement_ui <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "Chargement",
          h1("Subsample :"),
          
          #HTML("Additional text to disply only when menuItem tab One is selected \n"),
          #fileInput("dataFile",label = NULL, buttonLabel = "Parcourir...", placeholder = "Aucun patient sélectionné", accept = ".RData"),
          #selectInput("patient", "Sélectionnez le patient : ", choices = c("", "FL12C1888", "FL140304", "FL08G0293", "FL09C1164", "FULL"), selected = NULL, selectize = F),
          
          fluidRow(),
          
          fluidRow(
            column(width = 6,
                   box("Metadata filters : ", width = 12, solidHeader = T, collapsible = T,br(),br(),
                       column(6,wellPanel(radioButtons(inputId = "Subgroup", label = NULL, choices = metadata, selected = F))),
                       column(6,wellPanel(uiOutput("Dynamic_Sub_Spe")))
                   )),
            
            column(width = 3,
                   box("QC filters : ", width = 12, solidHeader = T, collapsible = T,
                       textInput("maximum", "Nombre maximal de gènes exprimé :", value = NULL, width = NULL, placeholder = NULL),
                       textInput("percent_mt", "Pourcentage seuil d'expression mitochondriale :", value = NULL, width = NULL, placeholder = NULL)
                   ),
                   
                   fluidRow(),
                   box("Gene expression filter : ", width = 12, solidHeader = T, collapsible = T,
                       selectInput('Svariables', "Gène : ", NULL, selectize=TRUE, selected = NULL),
                       textInput("Seuil_variables", "Expression supérieur à :", value = NULL, width = NULL, placeholder = NULL),
                   ),
                   
                   fluidRow(),
                   radioButtons(inputId = "combo", label = "Combinaison to analyse :", choiceNames = c("Excipient/RCHOP", "Pré-greffe/RCHOP", "Excipient/Pré-greffe"), choiceValues = c("EX_RCHOP","PG_RCHOP","PG_EX"),selected = F),
                   
                   fluidRow(),
                   div(actionButton(inputId = "actBtnPatient", label = "Subset",icon = icon("play") ), align = "left", style = "margin-bottom: 10px;", style = "margin-top: -10px;"),
                   div(actionButton(inputId = "resetPatient", label = "Reset",icon = icon("play") ),align = "left", style = "margin-bottom: 10px;", style = "margin-top: -10px;"),
            ))
)}

mod_Chargement_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_Chargement_ui("Chargement_ui_1")

## To be copied in the server
# callModule(mod_Chargement_server, "Chargement_ui_1")