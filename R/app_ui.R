#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd

options(shiny.autoload.r=FALSE)
app_ui <- function(request) {
  shinyUI(shinydashboard::dashboardPage(
    ######################### Titre ##########################
    shinydashboard::dashboardHeader(title = "Flinovo"),
    
    ########################## Menu ##########################
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        menuItem("Presentation", startExpanded = F,   tabName = "Presentation"),
        
        column(1,h4(textOutput("information"), style="color:white")),br(),br(),
        column(1,textOutput("nb_tot_cell")),br(),
        column(1,textOutput("nb_b_cell")),br(),
        column(1,textOutput("nb_other_cell")),br(),br(),
        column(1,textOutput("nb_Controle_cell")),br(),
        column(1,textOutput("nb_Excipient_cell")),br(),
        column(1,textOutput("nb_RCHOP_cell")),br(),br(),
        
        # All ------------------------------------------------
        menuItem("Meta-Analyses", startExpanded = T,
          
          column(1,h4("Dataset :", style="color:white")),br(),br(),
          box(width = 12,awesomeRadio(inputId = "Dataset_Type", label = "Type : ", choices = c("All","BT", "B", "Other"), selected = "All", inline = F, checkbox = TRUE)),br(),br(), #"BT"
          box(width = 12,awesomeRadio(inputId = "Dataset_Condition",label = "Condition : ",choices = c("All","Post-greffe","Pré-greffe-Excipient","Pré-greffe","Excipient","RCHOP"),selected = "All",inline = F,checkbox = TRUE)),br(),br(),br(),br(),br(),#br(),br(),br(),br(),br(),br(),br(),br(), #"All", "Pré-greffe", ,"RCHOP"
          actionButton(inputId = "actBtnLoadDataset", label = "Submit",icon = icon("play"), width ='87%'), br(),
          
          menuSubItem("Metadata",                     tabName = "All_Metadata", icon = icon("poll")),
          menuSubItem("Dimensional Reduction",        tabName = "All_Reduction_Dimension", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
          menuSubItem("Differentiel Expression",      tabName = "All_DE", icon = icon("calendar"))),
        #menuSubItem("Velocity",                     tabName = "All_Velocity", icon = icon("poll")))
        #menuSubItem("scRepertoir",                  tabName = "All_scRepertoir", icon = icon("poll")),
        #menuSubItem("Enrichissement",               tabName = "All_enrichissement", icon = icon("poll"))),

        
        # Patient -------------------------------------------
        menuItem("Analyses by patient", startExpanded = F,
                 selectInput("patient", h4("Patient selection : ", style="color:white"), choices = c("FL08G0293", "FL12C1888", "FL140304", "FL09C1164","FL02G095","FL05G0330","FL08G0431", "FL06G1206"), selected = "FL08G0293", selectize = T), #"FL1085","FL120316",  "FL1214",  "FL1481",
                 
                 box(width = 12,awesomeRadio(inputId = "Patient_Type", label = "Type : ", choices = c("All", "BT", "B"), selected = "All", inline = F, checkbox = TRUE)),br(),br(),
                 box(width = 12,awesomeRadio(inputId = "Patient_Condition",label = "Condition : ",choices = c("All","Pré-greffe-Excipient", "Post-greffe"),selected = "All",inline = F,checkbox = TRUE)), #"Pré-greffe"
                 actionButton(inputId = "actBtnLoadPatient", label = "Submit",icon = icon("play"), width ='87%'), br(),
                 
                 column(12,verbatimTextOutput("nb_clonotype_cell")),fluidRow(),br(),
                 
                 menuSubItem("Metadata",                     tabName = "Patient_Metadata", icon = icon("poll")),
                 menuSubItem("Dimensional Reduction",        tabName = "Patient_Reduction_Dimension", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
                 menuSubItem("VDJ Analyses",                 tabName = "Patient_VDJ", icon = icon("tasks", lib = "glyphicon")),
                 menuSubItem("ParTI Analyses",               tabName = "Patient_Parti", icon = icon("cog")),
                 menuSubItem("GoT Analyses",                 tabName = "Patient_GOT", icon = icon("knight", lib = "glyphicon")),
                 menuSubItem("Quality Control",              tabName = "Patient_QC", icon = icon("leaf", lib = "glyphicon")))
                #menuSubItem("Expression et enrichissement", tabName = "Patient_Expression", icon = icon("calendar"))),
                #menuSubItem("Subsamples",                   tabName = "Patient_Subset", icon = icon("cog"))),
      )
        
        

    ),
    
    ########################## Body ##########################
    shinydashboard::dashboardBody(
      ## -- Thème -- ##
      dashboardthemes::shinyDashboardThemes(theme = "blue_gradient"),
      tags$head(tags$style(HTML('.main-header .logo {font-family: "system-ui", Times, "Times New Roman", serif;font-weight: italic;font-size: 36px;'))),
      
      ## -- Page -- ##
      tabItems(
        mod_Presentation_ui("Presentation_ui_1"),
        
        mod_Patient_Metadata_ui("Patient_Metadata_ui_1"),
        mod_Patient_Reduction_Dimension_ui("Patient_Reduction_Dimension_ui_1"),
        mod_Patient_Parti_ui("Patient_Parti_ui_1"),
        mod_Patient_QC_ui("Patient_QC_ui_1"),
        mod_Patient_VDJ_ui("Patient_VDJ_ui_1"),
        mod_Patient_GOT_ui("Patient_GOT_ui_1"),
        #mod_Patient_Expression_ui("Patient_Expression_ui_1"),
        #mod_Patient_Subset_ui("Patient_Subset_ui_1"),
        
        
        mod_All_Metadata_ui("All_Metadata_ui_1"),
        mod_All_Reduction_Dimension_ui("All_Reduction_Dimension_ui_1"),
        mod_All_DE_ui("All_DE_ui_1")
        #mod_All_Cell_ui("All_Cell_ui_1"),
        #mod_All_Velocity_ui("All_Velocity_ui_1")
        #mod_All_scRepertoir_ui("All_scRepertoir_ui_1"),
        #mod_All_Enrichissement_ui("All_Enrichissement_ui_1"),
      ))
  ))
}
