mod_Patient_Metadata_ui <- function(id) {
  ns <- NS(id)
  # Métadonnées
  tabItem(tabName = "Patient_Metadata",
          h1("Présentation du patient et métadonnées : "),
          fluidPage(
            column(7,plotlyOutput(ns('dataTable'), width = "100%",  height = "800px")),
            column(5,
                verticalLayout(fluid = T,
                               tags$head(tags$style(HTML('.shiny-split-layout>div {overflow: hidden;}')),),
                               verbatimTextOutput(ns("nb_clonotype_cell2"))),
                downloadButton("PDF_report", "PDF report"),
                downloadButton("HTML_report", "HTML report"),
                downloadButton("Word_report", "Word report")
            )
          )
  )
}

mod_Patient_Metadata_server <- function(input, output, session, r) {
  ns <- session$ns
  
  
  ## -- Sunburst -- ##
  output$dataTable <- renderPlotly({
    plot_ly(r$dataset@tools$sunburst, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, marker = list(colors = c( "#BEBADA", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462")), type = 'sunburst',branchvalues = 'total', hoverinfo = "text", hovertext = paste(r$dataset@tools$sunburst$labels, ":", round((r$dataset@tools$sunburst$values/r$dataset@tools$sunburst$values[1])*100,2),"%", "\nTotal : " , r$dataset@tools$sunburst$values))
  })
  
  output$nb_clonotype_cell2 <- renderText({
    return(paste0("Numbre total de cellule détectée : \t", as.character(as.numeric(sum(table(r$dataset@meta.data$Phénotype)))),
                  "\n\nB cells : \t", as.character(as.numeric(table(r$dataset@meta.data$Phénotype)["B cells"])),
                  "\nOther cells : \t", as.character(as.numeric(sum(table(r$dataset@meta.data$Phénotype))) - as.numeric(table(r$dataset@meta.data$Phénotype)["B cells"])),
                  
                  "\n\nControle : \t", as.character(as.numeric(table(r$dataset@meta.data$Condition)["Pré-greffe"])),
                  "\nExcipient : \t", as.character(as.numeric(table(r$dataset@meta.data$Condition)["Excipient"])),
                  "\nRCHOP : \t", as.character(as.numeric(table(r$dataset@meta.data$Condition)["RCHOP"])),
                  
                  "\n\nClonotype majoritaire :", 
                  "\n Proportion : \t", round(r$dataset@tools$vloupe$proportion[1],3)*100,"% \n Type : \t", r$dataset@tools$vloupe$type[1],"\n Isotype : \t", r$dataset@tools$vloupe$igh_c_genes[1], 
                  "\n Heavy : \t", r$dataset@tools$vloupe$V_lourde[1], "/", r$dataset@tools$vloupe$D_lourde[1], "/", r$dataset@tools$vloupe$J_lourde[1], 
                  "\n Light : \t", r$dataset@tools$vloupe$V_legere[1], "/", r$dataset@tools$vloupe$J_legere[1], "\n"))
  })
  
  
}


## To be copied in the UI
# mod_Patient_Metadata_ui("Metadata_ui_1")

## To be copied in the server
# callModule(mod_Patient_Metadata_server, "Metadata_ui_1")