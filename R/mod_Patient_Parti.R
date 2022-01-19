mod_Patient_Parti_ui <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "Patient_Parti",
          fluidPage(
            fluidRow(
              column(9,
                conditionalPanel(condition = "input.Reduction == 'pca'", ns = ns,plotlyOutput(ns("Parti_PCA"), width = "100%",  height = "1000px")),
                conditionalPanel(condition = "input.Reduction == 'umap'", ns = ns,plotOutput(ns("Parti_UMAP"), width = "100%",  height = "1000px")),
                conditionalPanel(condition = "input.Reduction == 'tsne'", ns = ns,plotOutput(ns("Parti_TSNE"), width = "100%",  height = "1000px"))),
              column(3,
                shinyBS::bsCollapse(id = "collapse", open = "Paramètres : ",
                  shinyBS::bsCollapsePanel("Paramètres : ",
                     h5("Reduction :"),
                     pickerInput(inputId = ns("Reduction"), label = NULL , choices = c("pca","umap","tsne"), multiple = F, options = list(`actions-box` = TRUE), width = "100%"),fluidRow(),
                     h5("Archetype : "),
                     uiOutput(ns("archetypes")),fluidRow(),
                     h5("Condition :"),
                     pickerInput(inputId = ns("Condition"), label = NULL , choices = c("RCHOP","Excipient","Pré-greffe","Post-greffe","Pré-greffe-Excipient"), multiple = F, options = list(`actions-box` = TRUE), width = "100%"),fluidRow(),
                     h5("Database :"), 
                     pickerInput(inputId = ns("Database"), label = NULL , choices = c("CC","MF","BP"), multiple = F, options = list(`actions-box` = TRUE), width = "100%"),fluidRow(),
                     verbatimTextOutput(ns("archetype")),fluidRow(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br()
            ))))))
}

mod_Patient_Parti_server <- function(input, output, session, r, p) {
  ns <- session$ns
  reactive({input$Reduction})
  reactive({input$Condition})
  reactive({input$Archetype})
  reactive({input$Database})

  output$archetypes <- renderUI({
    load(file = paste0("datasets/ParTI/",input$Condition,"_",p$patient,"_parti_",input$Database,".RData"))
    pickerInput(inputId = ns("Archetype"),label = NULL , choices = c("Condition", "Phénotype", "Phénotype.fine", "Phase", "Clonotype","Chaine","V","D","J","C","CDR3", as.vector(labs$enriched_sets$y_name)),multiple = F,width = "100%", options = list(`actions-box` = TRUE))
  })
  output$archetype <- renderText({
    load(file = paste0("datasets/ParTI/",input$Condition,"_",p$patient,"_parti_",input$Database,".RData"))
    info <- stringr::str_replace(labs$lab$arch_lab, "\n\n", ": \n")
    info <- stringr::str_replace(info, "\n\n", "\n")
    info <- stringr::str_replace(info, "archetype_", "\narchetype ")
    info <- stringr::str_replace_all(info, "\n\n", "")
    info <- paste0(info,"\n")
    info[1] <- stringr::str_replace(info[1], "\n", "")
    info <- stringr::str_replace_all(info, "\n\n", "\n")
    return(info)
  })
  output$Parti_PCA <- renderPlotly({
    load(file = paste0("datasets/ParTI/",input$Condition,"_",p$patient,"_parti_",input$Database,".RData"))
    p_pca = ParetoTI::plot_arc(arc_data = arc, data = PCs4arch, which_dimensions = 1:3, data_lab = activ[[input$Archetype]], line_size = 3.5, text_size = 40, data_size = 3)
    plotly::layout(p_pca, title = paste(input$Archetype, "activity in", input$Condition))
  })
  output$Parti_UMAP <- renderPlot({
    load(file = paste0("datasets/ParTI/",input$Condition,"_",p$patient,"_parti_",input$Database,".RData"))
    arc_umap = ParetoTI::arch_to_umap(arc_data = arc, data = PCs4arch, which_dimensions = 1:2, method = c("naive", "umap-learn")) # implemented in R and slow, requires python module
    ParetoTI::plot_arc(arc_data = arc_umap$arc_data, data = arc_umap$data, which_dimensions = 1:2) + ggplot2::theme_bw()
  })
  output$Parti_TSNE <- renderPlot({
    load(file = paste0("datasets/ParTI/",input$Condition,"_",p$patient,"_parti_",input$Database,".RData"))
    arc_tsne = ParetoTI::arch_to_tsne(arc_data = arc, data = PCs4arch, which_dimensions = 1:2) # Project to tSNE coordinates (3D -> 2D, requires Rtsne package)
    ParetoTI::plot_arc(arc_data = arc_tsne$arc_data, data = arc_tsne$data, which_dimensions = 1:2) + ggplot2::theme_bw()
  })
}

## To be copied in the UI
# mod_Patient_Parti_ui("Patient_Parti_ui_1")

## To be copied in the server
# callModule(mod_Patient_Parti_server, "Patient_Parti_ui_1")