mod_Patient_VDJ_ui <- function(id) {
  ns <- NS(id)
  # VDJ
  tabItem(tabName = "Patient_VDJ",
          box("Clonotype", width = 12, plotlyOutput(ns("VDJ_Clonotype"), width = "100%",  height = "400px")),
          box("Heavy", width = 7, plotlyOutput(ns("VDJ_Heavy"), width = "100%",  height = "400px"),
              column(1,align="right",checkboxGroupInput("heavy_details", NULL, choices = list("Details" = "Details"), selected = 0)),
              conditionalPanel(
                condition = "input.heavy_details == 'Details' ", br(),
                plotlyOutput(ns("VDJ_DHeavy"), width = "100%",  height = "400px"))
          ),
          box("Lights", width = 5, plotlyOutput(ns("VDJ_Light"), width = "100%",  height = "400px"),
              column(1,align="right",checkboxGroupInput("light_details", NULL, choices = list("Details" = "Details"), selected = 0)),
              conditionalPanel(
                condition = "input.light_details == 'Details' ", br(),
                plotlyOutput(ns("VDJ_DLight"), width = "100%",  height = "400px"))
          ),
          fluidRow(
            box("V", width = 4, plotlyOutput(ns("V"), width = "100%",  height = "400px")),
            box("D", width = 4, plotlyOutput(ns("D"), width = "100%",  height = "400px")),
            box("J", width = 4, plotlyOutput(ns("J"), width = "100%",  height = "400px")))
  )
}

mod_Patient_VDJ_server <- function(input, output, session, r) {
  ns <- session$ns
  
  
  ## -- Clonotype -- ##
  output$VDJ_Clonotype = renderPlotly({
    plot_ly(x = r$dataset@tools$Clonotype[[1]], y = r$dataset@tools$Clonotype[[2]], name = "Info", type = "bar", 
            hovertemplate = paste0('Clonotype : %{x}\n', "Proportion : ", round(r$dataset@tools$vloupe$proportion[1:5],3)*100,"% \nType : ", r$dataset@tools$vloupe$type[1:5], "\nIsotype : ", r$dataset@tools$vloupe$igh_c_genes[1:5], 
                                   "\nHeavy : ", r$dataset@tools$vloupe$V_lourde[1:5], " / ", r$dataset@tools$vloupe$D_lourde[1:5], " / ", r$dataset@tools$vloupe$J_lourde[1:5], "\nLight : ", r$dataset@tools$vloupe$V_legere[1:5], " / ", r$dataset@tools$vloupe$J_legere[1:5])
    ) %>% layout(title='Frequencies of the mains clonotypes', yaxis =list(title="Number of cells"))
  })
  
  ## -- VDJ -- ##
  output$V = renderPlotly({plot_ly(x = r$dataset@tools$V[[1]], y = r$dataset@tools$V[[2]], name = "Clonotype", type = "bar") %>% layout(title='Frequencies V genes : Heavy and Lights chains', xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  output$D = renderPlotly({plot_ly(x = r$dataset@tools$D[[1]], y = r$dataset@tools$D[[2]], name = "Clonotype", type = "bar") %>% layout(title='Frequencies D genes : Heavy chain',xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  output$J = renderPlotly({plot_ly(x = r$dataset@tools$J[[1]], y = r$dataset@tools$J[[2]], name = "Clonotype", type = "bar") %>% layout(title='Frequencies J genes : Heavy and Lights chains',xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  
  ## -- Heavy chain -- ##
  output$VDJ_Heavy = renderPlotly({plot_ly(x = r$dataset@tools$Heavy[[1]], y = r$dataset@tools$Heavy[[2]], name = "Heavy Chain", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(r$dataset@tools$Heavy[[2]]/sum(r$dataset@tools$Heavy[[2]]),3)*100,"%")) %>% layout(title='Frequencies Heavy Chain' , xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))})
  output$VDJ_DHeavy = renderPlotly({plot_ly(x = r$dataset@tools$Isotype[[1]], y = r$dataset@tools$Isotype[[2]], name = "Heavy Chain", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(r$dataset@tools$Isotype[[2]]/sum(r$dataset@tools$Isotype[[2]]),3)*100,"%")) %>%layout(title='Frequencies Heavy Chain with details' , xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))  })
  
  ## -- Light chain -- ##
  output$VDJ_Light = renderPlotly({plot_ly(x = r$dataset@tools$Light[[1]], y = r$dataset@tools$Light[[2]], name = "Light Chain", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(r$dataset@tools$Light[[2]]/sum(r$dataset@tools$Light[[2]]),3)*100,"%")) %>% layout(title='Frequencies Light Chain',yaxis =list(title="Number of cells"))})
  output$VDJ_DLight = renderPlotly({plot_ly(x = r$dataset@tools$Type[[1]], y = r$dataset@tools$Type[[2]], name = "Clonotype", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(r$dataset@tools$Type[[2]]/sum(r$dataset@tools$Type[[2]]),3)*100,"%")) %>% layout(title='Frequencies Light Chain with details', yaxis =list(title="Number of cells"))})
  
  
}


## To be copied in the UI
# mod_Patient_VDJ_ui("VDJ_ui_1")

## To be copied in the server
# callModule(mod_Patient_VDJ_server, "VDJ_ui_1")