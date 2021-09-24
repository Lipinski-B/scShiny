mod_3D_ui <- function(id) {
  ns <- NS(id)
  # 3D Reduction dimention 
  tabItem(tabName = "3D_RD",
          radioButtons("metadata", "Metadata:",inline=T,singlet@tools$meta_variable),
          
          navbarPage("Reduction dimention",
                     tabPanel("PCA",plotlyOutput("D_PCA", width = "100%",  height = "800px"),uiOutput("Dfeature_pca")),
                     tabPanel("UMAP",plotlyOutput("D_UMAP", width = "100%",  height = "800px"),uiOutput("Dfeature_umap")),
                     tabPanel("TSNE",plotlyOutput("D_TSNE", width = "100%",  height = "800px"),uiOutput("Dfeature_tsne"))
          ),
          
          fluidRow(),
          column(1,align="right",checkboxGroupInput("DFeatures", NULL, choices = list("Features" = "Features"), selected = 0)),
          
          fluidRow(),
          conditionalPanel(
            condition = "input.DFeatures == 'Features'",
            fluidRow(column(6,selectInput('Dvariables', 'Choices', NULL, selectize=TRUE, selected = NULL))),
          ),
  )
}

mod_3D_server <- function(input, output, session) {
  ns <- session$ns
}


## To be copied in the UI
# mod_3D_ui("3D_ui_1")

## To be copied in the server
# callModule(mod_3D_server, "3D_ui_1")