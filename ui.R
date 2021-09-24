# ui.R
shinyUI(dashboardPage(
  ######################### Titre ##########################
  dashboardHeader(title = "Flinovo"),
  
  ########################## Menu ##########################
  dashboardSidebar(
    sidebarMenu(
      br(),
      selectInput("patient", h4("Sélectionnez le patient : ", style="color:white"), choices = c("FL08G0293", "FL12C1888", "FL140304", "FL09C1164","FL02G095","FL05G0330", "FL1085","FL120316",  "FL1214",  "FL1481","all"), selected = "FL08G0293", selectize = T),
      actionButton(inputId = "actBtnPatient1", label = "Submit",icon = icon("play"), width ='87%'),
      #submitButton("Submit", width = 230),
      br(),
      
      column(1,h4("Information :", style="color:white")),br(),br(), 
      column(1,textOutput("nb_tot_cell")),br(),
      column(1,textOutput("nb_b_cell")),br(),
      column(1,textOutput("nb_other_cell")),br(),br(), 
      column(1,textOutput("nb_Controle_cell")),br(),
      column(1,textOutput("nb_Excipient_cell")),br(),
      column(1,textOutput("nb_RCHOP_cell")),br(),br(),
      verbatimTextOutput("nb_clonotype_cell"),

      menuItem("Analyses", startExpanded = F,
        menuSubItem("Métadonnées",                  tabName = "metadata", icon = icon("poll")),
        menuSubItem("Réduction de dimension",       tabName = "readData", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
        #menuSubItem("3D Visualisation",            tabName = "3D_RD", icon = icon("indent-right", lib = "glyphicon")),
        #menuSubItem("Analyses PCA",                tabName = "PCA_facteurs", icon = icon("lock", lib = "glyphicon")),
        menuSubItem("Expression et enrichissement", tabName = "hallmark", icon = icon("calendar")),
        #menuSubItem("Trajectoires évolutives",     tabName = "monocle", icon = icon("sort-by-attributes-alt", lib = "glyphicon")),
        menuSubItem("Controle Qualité",             tabName = "mitochondrie", icon = icon("leaf", lib = "glyphicon")),
        menuSubItem("Analyses VDJ",                 tabName = "VDJ", icon = icon("readme")),
        menuSubItem("Subsamples",                   tabName = "Chargement", icon = icon("cog"))),
      menuItem("Méta-Analyses", startExpanded = F,
        menuSubItem("VDJ",                          tabName = "metadata_VDJ", icon = icon("poll")),
        menuSubItem("Cellules",                     tabName = "metadata_cell", icon = icon("poll")),
        menuSubItem("scRepertoir",                  tabName = "metadata_scRepertoir", icon = icon("poll")),
        menuSubItem("Enrichissement",               tabName = "metadata_enrichissement", icon = icon("poll")))
    )
  ),
  
  ########################## Body ##########################
  dashboardBody(
    ## -- Thème -- ##
    dashboardthemes::shinyDashboardThemes(theme = "blue_gradient"),
    tags$head(tags$style(HTML('.main-header .logo {font-family: "system-ui", Times, "Times New Roman", serif;font-weight: italic;font-size: 36px;'))),
    
    ## -- Page -- ##
    tabItems(# Chargement

      mod_Chargement_ui("Chargement_ui_1"),
      mod_Metadata_ui("Metadata_ui_1"),
      mod_Read_ui("Read_ui_1"),
      mod_PCA_ui("PCA_ui_1"),
      mod_Mitochondrie_ui("Mitochondrie_ui_1"),
      mod_Monocle_ui("Monocle_ui_1"),
      mod_Hallmark_ui("Hallmark_ui_1"),
      mod_VDJ_ui("VDJ_ui_1"),
      
      ############################################################################################
      
      mod_meta_VDJ_ui("meta_VDJ_ui_1"),
      mod_meta_Cell_ui("meta_Cell_ui_1"),
      mod_meta_scRepertoir_ui("meta_scRepertoir_ui_1"),
      mod_meta_Enrichissement_ui("meta_Enrichissement_ui_1")
      ))
))
